import numpy as np
import pandas as pd
from astropy.table import Table
import scipy.interpolate
from scipy.special import logsumexp

# Load spline track from Table 1
t1 = Table.read(r'D:\ARGON\desi_gd1_data\Jarvis_GD1_DESI_DR2_Zenodo\Table1.fits').to_pandas()
spl_dist = scipy.interpolate.CubicSpline(t1['phi1'], t1['distance'])

BASE = {
    'Q_t': 0.60, 'Q_c': 0.30,
    'mu_t': -0.29, 'mu_c': 0.72, 'sd_t': 2.49, 'sd_c': 6.13,
    'phi_t': -0.01, 'phi_c': 0.53, 'wphi_t': 0.23, 'wphi_c': 2.18,
    'pm1_t': 0.27, 'pm1_c': -5.78, 'wpm1_t': 2.65, 'wpm1_c': 41.36,
    'pm2_t': 0.39, 'pm2_c': -0.34, 'wpm2_t': 0.68, 'wpm2_c': 3.62,
}

def log_normal(x, mu, sd):
    return -0.5*np.log(2*np.pi) - np.log(sd) - 0.5*((x-mu)/sd)**2

cuts = [0.50, 0.60, 0.70, 0.80, 0.90, 0.95]

def setup_catalog(df):
    in_fit = (df['phi1'] > -60) & (df['phi1'] < 0)
    df_fit = df[in_fit].copy().reset_index(drop=True)
    N_stars = len(df_fit)
    
    d_arr = spl_dist(df_fit['phi1'].to_numpy(float))
    f = 4.74047 * d_arr
    f_ref = 4.74047 * 7.75
    
    xphi = df_fit['DELTA_PHI2'].to_numpy(float)
    xp1 = df_fit['DELTA_PM_PHI1'].to_numpy(float) * f
    xp2 = df_fit['DELTA_PM_PHI2'].to_numpy(float) * f
    ep1 = df_fit['PM_PHI1_ERR'].to_numpy(float) * f
    ep2 = df_fit['PM_PHI2_ERR'].to_numpy(float) * f
    ev = df_fit['V_ERR'].to_numpy(float)
    v_obs = df_fit['DELTA_VGSR'].to_numpy(float)
    
    ll_astro_t = (np.log(BASE['Q_t']) + log_normal(xphi, BASE['phi_t'], BASE['wphi_t'])
                  + log_normal(xp1, BASE['pm1_t'], np.hypot(BASE['wpm1_t'], ep1))
                  + log_normal(xp2, BASE['pm2_t'], np.hypot(BASE['wpm2_t'], ep2)))
    ll_astro_c = (np.log(BASE['Q_c']) + log_normal(xphi, BASE['phi_c'], BASE['wphi_c'])
                  + log_normal(xp1, BASE['pm1_c'], np.hypot(BASE['wpm1_c'], ep1))
                  + log_normal(xp2, BASE['pm2_c'], np.hypot(BASE['wpm2_c'], ep2)))
    ll_astro_b = np.full(N_stars, np.log(1 - BASE['Q_t'] - BASE['Q_c']) - np.log(12 * 4*f_ref * 4*f_ref * 60))
    
    sd_v_t = np.hypot(BASE['sd_t'], ev)
    sd_v_c = np.hypot(BASE['sd_c'], ev)
    
    posteriors = np.stack([df_fit['P_THIN'], df_fit['P_COCOON'], 1 - df_fit['P_THIN'] - df_fit['P_COCOON']], axis=1)
    posteriors = np.maximum(posteriors, 0)
    posteriors /= posteriors.sum(axis=1, keepdims=True)
    
    z_obs = np.abs((v_obs - BASE['mu_t']) / sd_v_t)
    
    return {
        'df': df_fit,
        'N_stars': N_stars,
        'll_astro_t': ll_astro_t,
        'll_astro_c': ll_astro_c,
        'll_astro_b': ll_astro_b,
        'sd_v_t': sd_v_t,
        'sd_v_c': sd_v_c,
        'posteriors': posteriors,
        'z_obs': z_obs
    }

def run_forward_sim(cat_dict, n_draws=10000, seed=42):
    rng = np.random.default_rng(seed)
    N_stars = cat_dict['N_stars']
    posteriors = cat_dict['posteriors']
    sd_v_t = cat_dict['sd_v_t']
    sd_v_c = cat_dict['sd_v_c']
    ll_astro_t = cat_dict['ll_astro_t']
    ll_astro_c = cat_dict['ll_astro_c']
    ll_astro_b = cat_dict['ll_astro_b']
    
    max_z_matrix = np.zeros((n_draws, len(cuts)))
    cum_p = np.cumsum(posteriors, axis=1)
    
    for b in range(n_draws):
        u = rng.uniform(0, 1, size=(N_stars, 1))
        comp = (u > cum_p[:, :2]).sum(axis=1)
        
        v_sim = np.zeros(N_stars)
        mask_t = (comp == 0)
        mask_c = (comp == 1)
        mask_b = (comp == 2)
        
        v_sim[mask_t] = rng.normal(BASE['mu_t'], sd_v_t[mask_t])
        v_sim[mask_c] = rng.normal(BASE['mu_c'], sd_v_c[mask_c])
        v_sim[mask_b] = rng.uniform(-30.0, 30.0, size=mask_b.sum())
        v_sim = np.clip(v_sim, -30.0, 30.0)
        
        ll_v_t = log_normal(v_sim, BASE['mu_t'], sd_v_t)
        ll_v_c = log_normal(v_sim, BASE['mu_c'], sd_v_c)
        
        LL = np.stack([ll_astro_t + ll_v_t, ll_astro_c + ll_v_c, ll_astro_b], axis=-1)
        P_sim = np.exp(LL - logsumexp(LL, axis=-1)[:, None])
        p_thin_sim = P_sim[:, 0]
        z_sim = np.abs((v_sim - BASE['mu_t']) / sd_v_t)
        
        for j, c in enumerate(cuts):
            pass_c = p_thin_sim >= c
            if pass_c.sum() > 0:
                max_z_matrix[b, j] = np.max(z_sim[pass_c])
                
    return max_z_matrix

def compute_look_elsewhere(M, obs_max_z):
    n_draws, n_cuts = M.shape
    p_pointwise = np.mean(M >= obs_max_z[None, :], axis=0)
    p_min_obs = np.min(p_pointwise)
    
    sim_p_matrix = np.zeros_like(M)
    for j in range(n_cuts):
        col = M[:, j]
        # empirical sf
        sim_p_matrix[:, j] = np.mean(col[:, None] <= col[None, :], axis=1)
        
    sim_min_p = np.min(sim_p_matrix, axis=1)
    p_global = np.mean(sim_min_p <= p_min_obs)
    
    return p_pointwise, p_min_obs, p_global

print("Loading Zenodo v3 and v4...")
df_v3 = pd.read_csv('D:/ARGON/Table7.csv', dtype={'SOURCE_ID': str})
df_v4 = Table.read('D:/ARGON/desi_gd1_data/v4_extracted/Jarvis_GD1_DESI_DR2_Zenodo/Table7_all.fits').to_pandas()
df_v4['SOURCE_ID'] = df_v4['SOURCE_ID'].astype(str).str.strip()

cat_v3 = setup_catalog(df_v3)
cat_v4 = setup_catalog(df_v4)

print("Running 10,000 draws on v3...")
M_v3 = run_forward_sim(cat_v3, n_draws=10000, seed=42)

print("Running 10,000 draws on v4...")
M_v4 = run_forward_sim(cat_v4, n_draws=10000, seed=42)

# V3 observed
s28_idx_v3 = np.where(cat_v3['df']['SOURCE_ID'] == '804436551223805568')[0][0]
obs_max_z_v3_full = np.array([np.max(cat_v3['z_obs'][cat_v3['df']['P_THIN'] >= c]) for c in cuts])
obs_max_z_v3_no28 = np.array([np.max(cat_v3['z_obs'][(cat_v3['df']['P_THIN'] >= c) & (cat_v3['df'].index != s28_idx_v3)]) for c in cuts])

pt_v3_full, min_p_v3_full, glob_v3_full = compute_look_elsewhere(M_v3, obs_max_z_v3_full)
pt_v3_no28, min_p_v3_no28, glob_v3_no28 = compute_look_elsewhere(M_v3, obs_max_z_v3_no28)

# V4 observed
s28_idx_v4 = np.where(cat_v4['df']['SOURCE_ID'] == '804436551223805568')[0][0]
obs_max_z_v4_full = np.array([np.max(cat_v4['z_obs'][cat_v4['df']['P_THIN'] >= c]) for c in cuts])
obs_max_z_v4_no28 = np.array([np.max(cat_v4['z_obs'][(cat_v4['df']['P_THIN'] >= c) & (cat_v4['df'].index != s28_idx_v4)]) for c in cuts])

pt_v4_full, min_p_v4_full, glob_v4_full = compute_look_elsewhere(M_v4, obs_max_z_v4_full)
pt_v4_no28, min_p_v4_no28, glob_v4_no28 = compute_look_elsewhere(M_v4, obs_max_z_v4_no28)

print("\n==================== ZENODO V3 RESULTS ====================")
print("--- Full Sample (v3) ---")
for c, obs, p_pt in zip(cuts, obs_max_z_v3_full, pt_v3_full):
    print(f"P >= {c:.2f}: Obs Max |z| = {obs:.4f}, Pointwise p = {p_pt:.4f}")
print(f"Minimum Pointwise p: {min_p_v3_full:.6f}")
print(f"Global p-value (Look-Elsewhere Corrected): {glob_v3_full:.4f}")

print("\n--- Excluding Star 28 (v3) ---")
for c, obs, p_pt in zip(cuts, obs_max_z_v3_no28, pt_v3_no28):
    print(f"P >= {c:.2f}: Obs Max |z| = {obs:.4f}, Pointwise p = {p_pt:.4f}")
print(f"Minimum Pointwise p (no 28): {min_p_v3_no28:.6f}")
print(f"Global p-value (Look-Elsewhere Corrected): {glob_v3_no28:.4f}")

print("\n==================== ZENODO V4 RESULTS ====================")
print("--- Full Sample (v4) ---")
for c, obs, p_pt in zip(cuts, obs_max_z_v4_full, pt_v4_full):
    print(f"P >= {c:.2f}: Obs Max |z| = {obs:.4f}, Pointwise p = {p_pt:.4f}")
print(f"Minimum Pointwise p: {min_p_v4_full:.6f}")
print(f"Global p-value (Look-Elsewhere Corrected): {glob_v4_full:.4f}")

print("\n--- Excluding Star 28 (v4) ---")
for c, obs, p_pt in zip(cuts, obs_max_z_v4_no28, pt_v4_no28):
    print(f"P >= {c:.2f}: Obs Max |z| = {obs:.4f}, Pointwise p = {p_pt:.4f}")
print(f"Minimum Pointwise p (no 28): {min_p_v4_no28:.6f}")
print(f"Global p-value (Look-Elsewhere Corrected): {glob_v4_no28:.4f}")
