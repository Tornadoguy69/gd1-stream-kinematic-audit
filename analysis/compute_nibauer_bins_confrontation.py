import numpy as np
import pandas as pd
from astropy.table import Table
from scipy.optimize import minimize_scalar

print("=" * 80)
print("BIN-BY-BIN COMPARISON WITH NIBAUER ET AL. (2026) BINS")
print("Bins: [-60, -40], [-40, -20], [-20, 0]")
print("=" * 80)

# Load Table7_all.fits
t = Table.read('D:/ARGON/Table7_all.fits').to_pandas()
t['SOURCE_ID'] = t['SOURCE_ID'].astype(str).str.strip()

# 9 candidate binary IDs from v45
type_a_ids = {
    '804436551223805568', # Star 28 (phi1 = -32.93)
    '637140290779827072', # phi1 = -57.93
    '1577287126973726592', # phi1 = -3.25
    '794510439424008832', # phi1 = -42.26
    '829645264771612800', # Star 4 (phi1 = -27.35)
    '794347947923184640', # phi1 = -43.79
    '839367829479258240', # phi1 = -18.48
    '1577024756011421696', # phi1 = -2.74
    '796893149840353408', # Outlier 13 (phi1 = -39.46)
}

# Profile likelihood function for intrinsic dispersion in a bin
def fit_sigma_int(df, weight_col='P_THIN'):
    w = df[weight_col].values
    v = df['DELTA_VGSR'].values
    ev = df['V_ERR'].values
    
    # Grid search or minimize negative log-likelihood
    def nll(sig):
        # Center mu
        var = sig**2 + ev**2
        mu = np.sum(w * v / var) / np.sum(w / var)
        term1 = np.log(2 * np.pi * var)
        term2 = (v - mu)**2 / var
        return 0.5 * np.sum(w * (term1 + term2))
    
    res = minimize_scalar(nll, bounds=(0.1, 10.0), method='bounded')
    best_sig = res.x
    
    # Standard bootstrap error
    boot_sigs = []
    n = len(df)
    for _ in range(500):
        idx = np.random.choice(n, n, replace=True)
        w_b, v_b, ev_b = w[idx], v[idx], ev[idx]
        def nll_b(sig):
            var = sig**2 + ev_b**2
            mu = np.sum(w_b * v_b / var) / np.sum(w_b / var)
            return 0.5 * np.sum(w_b * (np.log(2 * np.pi * var) + (v_b - mu)**2 / var))
        r_b = minimize_scalar(nll_b, bounds=(0.1, 10.0), method='bounded')
        boot_sigs.append(r_b.x)
    
    sig_err = np.std(boot_sigs)
    return best_sig, sig_err

bins = [
    ("Trailing Bin", -60.0, -40.0),
    ("Hottest / Spur Bin", -40.0, -20.0),
    ("Leading Bin", -20.0, 0.0)
]

results = []

for name, p_min, p_max in bins:
    bin_all = t[(t['phi1'] >= p_min) & (t['phi1'] < p_max) & ((t['P_THIN'] + t['P_COCOON']) > 0.5)].copy()
    bin_clean = bin_all[~bin_all['SOURCE_ID'].isin(type_a_ids)].copy()
    
    sig_all, err_all = fit_sigma_int(bin_all)
    sig_clean, err_clean = fit_sigma_int(bin_clean)
    
    n_a_in_bin = bin_all['SOURCE_ID'].isin(type_a_ids).sum()
    
    results.append({
        'bin': name,
        'phi1_range': f"[{p_min:.0f}, {p_max:.0f}]",
        'N_raw': len(bin_all),
        'N_binaries': int(n_a_in_bin),
        'N_clean': len(bin_clean),
        'sigma_raw': sig_all,
        'sigma_raw_err': err_all,
        'sigma_clean': sig_clean,
        'sigma_clean_err': err_clean,
        'reduction_pct': (sig_all - sig_clean) / sig_all * 100.0
    })

res_df = pd.DataFrame(results)
print(res_df.to_string(index=False))
res_df.to_csv("D:/ARGON/nibauer_bin_comparison.csv", index=False)
