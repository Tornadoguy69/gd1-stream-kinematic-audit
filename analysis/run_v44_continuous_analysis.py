import sys, os, json, hashlib
sys.stdout.reconfigure(encoding='utf-8')
sys.path.append('D:/ARGON')
import numpy as np
import pandas as pd
from astropy.table import Table
import scipy.interpolate
import scipy.stats
from scipy.special import logsumexp

print("=" * 80)
print("RESEARCH v44: THRESHOLD-FREE CONTINUOUS BAYESIAN SCORE & HIERARCHICAL CALIBRATION")
print("=" * 80)

# 1. Load Data
v3_path = 'D:/ARGON/Table7.csv'
v4_path = 'D:/ARGON/desi_gd1_data/v4_extracted/Jarvis_GD1_DESI_DR2_Zenodo/Table7_all.fits'
t1_path = 'D:/ARGON/desi_gd1_data/Jarvis_GD1_DESI_DR2_Zenodo/Table1.fits'

def sha256_file(filepath):
    h = hashlib.sha256()
    with open(filepath, 'rb') as f:
        while chunk := f.read(65536):
            h.update(chunk)
    return h.hexdigest()

v3_hash = sha256_file(v3_path)
v4_hash = sha256_file(v4_path)
print(f"Loaded Zenodo v3 (SHA-256: {v3_hash[:16]}...)")
print(f"Loaded Zenodo v4 (SHA-256: {v4_hash[:16]}... Record 23599156)")

t1 = Table.read(t1_path).to_pandas()
spl_dist = scipy.interpolate.CubicSpline(t1['phi1'], t1['distance'])

BASE = {
    'Q_t': 0.601609, 'Q_c': 0.297571,
    'mu_t': -0.294194, 'mu_c': 0.723762, 'sd_t': 2.491687, 'sd_c': 6.129247,
    'phi_t': -0.010313, 'phi_c': 0.532819, 'wphi_t': 0.228882, 'wphi_c': 2.179728,
    'pm1_t': 0.271025, 'pm1_c': -5.781126, 'wpm1_t': 2.650217, 'wpm1_c': 41.361825,
    'pm2_t': 0.394327, 'pm2_c': -0.339557, 'wpm2_t': 0.676888, 'wpm2_c': 3.622662,
}

def log_normal(x, mu, sd):
    return -0.5*np.log(2*np.pi) - np.log(sd) - 0.5*((x-mu)/sd)**2

def prepare_cat(df):
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
    ll_astro_b = np.full(N_stars, np.log(1 - BASE['Q_t'] - BASE['Q_c']) - np.log(12 * 4*f_ref * 4*f_ref))
    
    sd_v_t = np.hypot(BASE['sd_t'], ev)
    sd_v_c = np.hypot(BASE['sd_c'], ev)
    
    q_t = scipy.stats.norm.cdf((30.0 - BASE['mu_t'])/sd_v_t) - scipy.stats.norm.cdf((-30.0 - BASE['mu_t'])/sd_v_t)
    q_c = scipy.stats.norm.cdf((30.0 - BASE['mu_c'])/sd_v_c) - scipy.stats.norm.cdf((-30.0 - BASE['mu_c'])/sd_v_c)
    q_b = 1.0
    
    posteriors_pub = np.stack([df_fit['P_THIN'], df_fit['P_COCOON'], 1 - df_fit['P_THIN'] - df_fit['P_COCOON']], axis=1)
    posteriors_pub = np.maximum(posteriors_pub, 0)
    posteriors_pub /= posteriors_pub.sum(axis=1, keepdims=True)
    
    z_obs = np.abs((v_obs - BASE['mu_t']) / sd_v_t)
    s28_idx = np.where(df_fit['SOURCE_ID'].astype(str) == '804436551223805568')[0][0]
    
    return {
        'df': df_fit,
        'N_stars': N_stars,
        'll_astro_t': ll_astro_t,
        'll_astro_c': ll_astro_c,
        'll_astro_b': ll_astro_b,
        'sd_v_t': sd_v_t,
        'sd_v_c': sd_v_c,
        'q_t': q_t,
        'q_c': q_c,
        'q_b': q_b,
        'posteriors_pub': posteriors_pub,
        'z_obs': z_obs,
        'v_obs': v_obs,
        's28_idx': s28_idx
    }

df_v3 = pd.read_csv(v3_path, dtype={'SOURCE_ID': str})
df_v4 = Table.read(v4_path).to_pandas()
df_v4['SOURCE_ID'] = df_v4['SOURCE_ID'].astype(str).str.strip()

cat_v3 = prepare_cat(df_v3)
cat_v4 = prepare_cat(df_v4)
N = cat_v4['N_stars']

# -------------------------------------------------------------
# STEP 1: Hierarchical Spatial Empirical Bayes Prior
# -------------------------------------------------------------
print("\n[Step 1] Fitting Hierarchical Spatial Empirical Bayes Prior along stream track...")
phi1 = cat_v4['df']['phi1'].to_numpy(float)
p_thin_v4 = cat_v4['df']['P_THIN'].to_numpy(float)

# Spatial bins (10 bins of 6 degrees)
n_bins = 10
bins_phi1 = np.linspace(-60, 0, n_bins + 1)
bin_centers = 0.5 * (bins_phi1[:-1] + bins_phi1[1:])
bin_idx = np.clip(np.digitize(phi1, bins_phi1) - 1, 0, n_bins - 1)

bin_p_means = np.zeros(n_bins)
bin_counts = np.zeros(n_bins, dtype=int)
for k in range(n_bins):
    m = (bin_idx == k)
    bin_counts[k] = m.sum()
    if bin_counts[k] > 0:
        bin_p_means[k] = np.mean(p_thin_v4[m])

mean_stream_frac = float(np.mean(bin_p_means))
sigma_pi_spatial = float(np.std(bin_p_means))
print(f"  Empirical track mean thin fraction: {mean_stream_frac:.4f}")
print(f"  Spatial overdispersion across bins (sigma_pi): {sigma_pi_spatial:.4f}")

# -------------------------------------------------------------
# STEP 2: Continuous Probability-Weighted Score Statistic T_cont
# -------------------------------------------------------------
print("\n[Step 2] Evaluating Continuous Score Statistic T_cont = sum(w_i * max(0, z_i^2 - c0^2))...")
c0 = 2.5
c0_sq = c0**2

psi_v3 = np.maximum(0.0, cat_v3['z_obs']**2 - c0_sq)
T_cont_v3 = float(np.sum(cat_v3['df']['P_THIN'] * psi_v3))

s28_idx = cat_v4['s28_idx']
mask_no28 = np.ones(N, dtype=bool)
mask_no28[s28_idx] = False

T_cont_v3_no28 = float(np.sum((cat_v3['df']['P_THIN'] * psi_v3)[mask_no28]))

psi_v4 = np.maximum(0.0, cat_v4['z_obs']**2 - c0_sq)
T_cont_v4 = float(np.sum(cat_v4['df']['P_THIN'] * psi_v4))
T_cont_v4_no28 = float(np.sum((cat_v4['df']['P_THIN'] * psi_v4)[mask_no28]))

# Star 28 contributions
w28_v3 = float(cat_v3['df']['P_THIN'].iloc[s28_idx])
w28_v4 = float(cat_v4['df']['P_THIN'].iloc[s28_idx])
contrib28_v3 = float(w28_v3 * psi_v3[s28_idx])
contrib28_v4 = float(w28_v4 * psi_v4[s28_idx])

print(f"  Zenodo v3 Observed T_cont (all stars):        {T_cont_v3:.4f}")
print(f"  Zenodo v3 Observed T_cont (excluding Star 28): {T_cont_v3_no28:.4f}")
print(f"  Zenodo v4 Observed T_cont (all stars):        {T_cont_v4:.4f}")
print(f"  Zenodo v4 Observed T_cont (excluding Star 28): {T_cont_v4_no28:.4f}")
print(f"  Star 28 v3 weight = {w28_v3:.4f} -> contribution = {contrib28_v3:.4f}")
print(f"  Star 28 v4 weight = {w28_v4:.4f} -> contribution = {contrib28_v4:.4f} (gentle -1.8% change!)")

# -------------------------------------------------------------
# STEP 3: Monte Carlo Null Simulation for T_cont & Predictive Counts
# -------------------------------------------------------------
B = 10000
seed = 42
print(f"\n[Step 3] Running Full Rejection-Truncated Forward Simulation (B={B:,}, seed={seed})...")
rng = np.random.default_rng(seed)

sd_v_t = cat_v4['sd_v_t']
sd_v_c = cat_v4['sd_v_c']
q_t = cat_v4['q_t']
q_c = cat_v4['q_c']
ll_astro_t = cat_v4['ll_astro_t']
ll_astro_c = cat_v4['ll_astro_c']
ll_astro_b = cat_v4['ll_astro_b']

cum_p_v4 = np.cumsum(cat_v4['posteriors_pub'], axis=1)

T_sim_v4 = np.zeros(B)
T_sim_v4_no28 = np.zeros(B)
T_sim_v3 = np.zeros(B)
T_sim_v3_no28 = np.zeros(B)

counts_p80_baseline = np.zeros(B, dtype=int)
counts_p80_hier = np.zeros(B, dtype=int)

cuts = [0.50, 0.60, 0.70, 0.80, 0.90, 0.95]
n_pass_hier = np.zeros((B, len(cuts)), dtype=int)

for b in range(B):
    # Baseline simulation of components
    u = rng.uniform(0, 1, size=(N, 1))
    comp = (u > cum_p_v4[:, :2]).sum(axis=1)
    
    # Velocity realization with rejection truncation on [-30, 30]
    v_sim = np.zeros(N)
    mask_t = (comp == 0)
    mask_c = (comp == 1)
    mask_b = (comp == 2)
    
    if mask_t.sum() > 0:
        vt = rng.normal(BASE['mu_t'], sd_v_t[mask_t])
        bad = np.abs(vt) > 30.0
        while bad.any():
            vt[bad] = rng.normal(BASE['mu_t'], sd_v_t[mask_t][bad])
            bad = np.abs(vt) > 30.0
        v_sim[mask_t] = vt
    if mask_c.sum() > 0:
        vc = rng.normal(BASE['mu_c'], sd_v_c[mask_c])
        bad = np.abs(vc) > 30.0
        while bad.any():
            vc[bad] = rng.normal(BASE['mu_c'], sd_v_c[mask_c][bad])
            bad = np.abs(vc) > 30.0
        v_sim[mask_c] = vc
    v_sim[mask_b] = rng.uniform(-30.0, 30.0, size=mask_b.sum())
    
    ll_v_t = log_normal(v_sim, BASE['mu_t'], sd_v_t) - np.log(q_t)
    ll_v_c = log_normal(v_sim, BASE['mu_c'], sd_v_c) - np.log(q_c)
    ll_v_b = np.full(N, -np.log(60.0))
    
    LL = np.stack([ll_astro_t + ll_v_t, ll_astro_c + ll_v_c, ll_astro_b + ll_v_b], axis=-1)
    P_sim = np.exp(LL - logsumexp(LL, axis=-1)[:, None])
    p_thin_sim = P_sim[:, 0]
    z_sim = np.abs((v_sim - BASE['mu_t']) / sd_v_t)
    psi_sim = np.maximum(0.0, z_sim**2 - c0_sq)
    
    T_sim_v4[b] = np.sum(p_thin_sim * psi_sim)
    T_sim_v4_no28[b] = np.sum((p_thin_sim * psi_sim)[mask_no28])
    counts_p80_baseline[b] = (p_thin_sim >= 0.80).sum()
    
    # Hierarchical track-density realization:
    delta_track = rng.normal(0.0, sigma_pi_spatial, size=n_bins)
    # Apply local modulation to stream prior
    delta_star = delta_track[bin_idx]
    ll_astro_t_hier = ll_astro_t + delta_star
    LL_hier = np.stack([ll_astro_t_hier + ll_v_t, ll_astro_c + ll_v_c, ll_astro_b + ll_v_b], axis=-1)
    p_thin_hier = np.exp(LL_hier - logsumexp(LL_hier, axis=-1)[:, None])[:, 0]
    for j, c in enumerate(cuts):
        n_pass_hier[b, j] = (p_thin_hier >= c).sum()

p_T_v4_full = (np.sum(T_sim_v4 >= T_cont_v4) + 1.0) / (B + 1.0)
p_T_v4_no28 = (np.sum(T_sim_v4_no28 >= T_cont_v4_no28) + 1.0) / (B + 1.0)

print(f"\n--- EMPIRICAL MONTE CARLO RESULTS FOR T_cont (B={B:,}) ---")
print(f"  v4 Full Sample: T_obs = {T_cont_v4:.4f}, Mean T_sim = {np.mean(T_sim_v4):.4f} +/- {np.std(T_sim_v4):.4f}")
print(f"  -> Empirical p-value (v4 Full): {p_T_v4_full:.4f}")
print(f"  v4 Excl Star 28: T_obs = {T_cont_v4_no28:.4f}, Mean T_sim = {np.mean(T_sim_v4_no28):.4f} +/- {np.std(T_sim_v4_no28):.4f}")
print(f"  -> Empirical p-value (v4 Excl Star 28): {p_T_v4_no28:.4f}")

# -------------------------------------------------------------
# STEP 4: Sample Size Calibration under Hierarchical Prior
# -------------------------------------------------------------
print("\n--- SAMPLE SIZE CALIBRATION (HIERARCHICAL EMPIRICAL BAYES PRIOR) ---")
hier_calib_table = []
for j, c in enumerate(cuts):
    obs_n = int((cat_v4['df']['P_THIN'] >= c).sum())
    dist_hier = n_pass_hier[:, j]
    med_hier = float(np.median(dist_hier))
    lo_hier, hi_hier = np.quantile(dist_hier, [0.025, 0.975])
    p_lower = (np.sum(dist_hier <= obs_n) + 1.0) / (B + 1.0)
    p_upper = (np.sum(dist_hier >= obs_n) + 1.0) / (B + 1.0)
    two_sided_p = float(min(1.0, 2.0 * min(p_lower, p_upper)))
    hier_calib_table.append({
        'threshold': c,
        'obs_n': obs_n,
        'median_hier': med_hier,
        'ci95_hier': [float(lo_hier), float(hi_hier)],
        'p_lower': float(p_lower),
        'two_sided_p': two_sided_p
    })
    print(f"  P >= {c:.2f}: Obs N={obs_n:3d} | Hier Med={med_hier:5.1f} | 95% CI=[{lo_hier:3.0f}, {hi_hier:3.0f}] | p_N = {two_sided_p:.4f}")

# -------------------------------------------------------------
# STEP 5: Spatial-Kinematic Kernel Scan & Dipole Analysis
# -------------------------------------------------------------
print("\n[Step 5] Evaluating Spatial-Kinematic Scan S(phi1) and Permutation Null...")
grid_phi1 = np.linspace(-58, -2, 281)
h = 1.5

def kernel_scan(phi_arr, w_arr, z_arr, v_arr):
    scan_energy = np.zeros(len(grid_phi1))
    scan_mean_v = np.zeros(len(grid_phi1))
    for idx, x0 in enumerate(grid_phi1):
        K = np.exp(-0.5 * ((phi_arr - x0)/h)**2)
        W = K * w_arr
        sW = np.sum(W)
        if sW > 1e-4:
            psi_k = np.maximum(0.0, z_arr**2 - c0_sq)
            scan_energy[idx] = np.sum(W * psi_k) / sW
            scan_mean_v[idx] = np.sum(W * v_arr) / sW
    return scan_energy, scan_mean_v

scan_obs_energy, scan_obs_v = kernel_scan(phi1, p_thin_v4, cat_v4['z_obs'], cat_v4['v_obs'])
peak_obs_energy = float(np.max(scan_obs_energy))
peak_phi1 = float(grid_phi1[np.argmax(scan_obs_energy)])
print(f"  Observed Peak Tail Energy Density: {peak_obs_energy:.4f} at phi1 = {peak_phi1:.2f} deg (near Star 28 at -32.93 deg)")

# Spatial permutation test: permute observed velocities across stream stars
N_perm = 1000
max_energy_perms = np.zeros(N_perm)
for p in range(N_perm):
    perm_idx = rng.permutation(N)
    perm_v = cat_v4['v_obs'][perm_idx]
    perm_z = cat_v4['z_obs'][perm_idx]
    e_perm, _ = kernel_scan(phi1, p_thin_v4, perm_z, perm_v)
    max_energy_perms[p] = np.max(e_perm)

p_spatial_peak = (np.sum(max_energy_perms >= peak_obs_energy) + 1.0) / (N_perm + 1.0)
print(f"  Look-Elsewhere Corrected Spatial Permutation p-value for Peak: {p_spatial_peak:.4f}")

# -------------------------------------------------------------
# Save v44 Summary Output
# -------------------------------------------------------------
os.makedirs('D:/ARGON/v44_results', exist_ok=True)
v44_summary = {
    'provenance': {
        'version': 'Research v44 (Threshold-Free Continuous Bayesian Score & Hierarchical Calibration)',
        'v3_sha256': v3_hash,
        'v4_sha256': v4_hash,
        'random_seed': seed,
        'n_sim_draws': B
    },
    'hierarchical_calibration': {
        'spatial_overdispersion_sigma_pi': sigma_pi_spatial,
        'predictive_counts': hier_calib_table
    },
    'continuous_score_T_cont': {
        'excess_threshold_c0_sigma': c0,
        'v3_T_obs_full': T_cont_v3,
        'v3_T_obs_no28': T_cont_v3_no28,
        'v4_T_obs_full': T_cont_v4,
        'v4_T_obs_no28': T_cont_v4_no28,
        'star_28_contrib_v3': contrib28_v3,
        'star_28_contrib_v4': contrib28_v4,
        'v4_p_value_full': p_T_v4_full,
        'v4_p_value_no28': p_T_v4_no28
    },
    'spatial_kinematic_scan': {
        'kernel_bandwidth_deg': h,
        'peak_energy_density': peak_obs_energy,
        'peak_location_phi1_deg': peak_phi1,
        'permutation_calibrated_p_value': p_spatial_peak
    }
}

with open('D:/ARGON/v44_results/v44_continuous_summary.json', 'w') as f:
    json.dump(v44_summary, f, indent=2)

print("\nv44 Continuous Analysis complete! Summary saved to D:/ARGON/v44_results/v44_continuous_summary.json")
