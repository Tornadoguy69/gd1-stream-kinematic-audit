import urllib.request
import urllib.parse
import json
import numpy as np
import pandas as pd
from astropy.table import Table

print("=" * 80)
print("EMPIRICAL SDSS DR17 vs DESI DR2 ZERO-POINT OFFSET MEASUREMENT")
print("=" * 80)

# Load DESI stars
t = Table.read('D:/ARGON/Table7_all.fits').to_pandas()
sub = t[(t['phi1'] > -60) & (t['phi1'] < 0)].reset_index(drop=True)
print(f"Total DESI stars in fit window: {len(sub)}")

batch_size = 20
matches = []

for b_start in range(0, len(sub), batch_size):
    batch = sub.iloc[b_start : b_start + batch_size]
    conds = []
    for _, row in batch.iterrows():
        ra = row['RA']
        dec = row['Dec']
        r = 2.0 / 3600.0 # 2 arcseconds
        cosd = np.cos(np.radians(dec))
        conds.append(f"(ra BETWEEN {ra - r/cosd:.6f} AND {ra + r/cosd:.6f} AND dec BETWEEN {dec - r:.6f} AND {dec + r:.6f})")
    
    clause = " OR ".join(conds)
    sql = f"""
    SELECT specObjID, bestObjID, ra, dec, mjd, plate, fiberID, zwarning,
           (z * 299792.458) as v_sdss, (zErr * 299792.458) as ev_sdss
    FROM specObjAll
    WHERE zwarning = 0 AND ({clause})
    """
    url = "https://skyserver.sdss.org/dr17/SkyServerWS/SearchTools/SqlSearch?cmd=" + urllib.parse.quote(sql) + "&format=json"
    req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64)'})
    
    try:
        with urllib.request.urlopen(req, timeout=15) as resp:
            data = json.loads(resp.read().decode('utf-8'))
            rows = data[0].get('Rows', [])
            for r in rows:
                ra_s = float(r['ra'])
                dec_s = float(r['dec'])
                # match to closest DESI star
                cosd_s = np.cos(np.radians(dec_s))
                dists = np.hypot((batch['RA'] - ra_s) * cosd_s, batch['Dec'] - dec_s) * 3600.0
                closest_idx = dists.idxmin()
                if dists[closest_idx] <= 2.5: # within 2.5 arcsec
                    d_row = batch.loc[closest_idx]
                    v_desi = float(d_row['V_LOS'])
                    ev_desi = float(d_row['V_ERR'])
                    v_sdss = float(r['v_sdss'])
                    ev_sdss = float(r['ev_sdss'])
                    matches.append({
                        'source_id': str(d_row['SOURCE_ID']),
                        'phi1': float(d_row['phi1']),
                        'P_thin': float(d_row['P_THIN']),
                        'v_desi': v_desi,
                        'ev_desi': ev_desi,
                        'v_sdss': v_sdss,
                        'ev_sdss': ev_sdss,
                        'delta_v': v_desi - v_sdss,
                        'sigma_diff': np.hypot(ev_desi, ev_sdss),
                        'specObjID': str(r['specObjID']),
                        'plate': int(r['plate']),
                        'fiber': int(r['fiberID']),
                        'mjd': float(r['mjd']),
                        'dist_arcsec': float(dists[closest_idx])
                    })
    except Exception as e:
        print(f"Batch {b_start//batch_size} failed: {e}")

m_df = pd.DataFrame(matches)
# Drop duplicate specObjID matches if any
m_df = m_df.drop_duplicates(subset=['source_id', 'specObjID']).reset_index(drop=True)

print(f"\nTotal empirical cross-matches found: {len(m_df)}")
m_df.to_csv("D:/ARGON/sdss_desi_empirical_crossmatches.csv", index=False)

if len(m_df) > 0:
    # All stars
    med_all = np.median(m_df['delta_v'])
    mean_all = np.mean(m_df['delta_v'])
    std_all = np.std(m_df['delta_v'])
    
    # Exclude outliers (|delta_v| > 15 km/s, e.g. binaries)
    clean = m_df[m_df['delta_v'].abs() < 15.0]
    med_clean = np.median(clean['delta_v'])
    mean_clean = np.mean(clean['delta_v'])
    std_clean = np.std(clean['delta_v'])
    sem_clean = std_clean / np.sqrt(len(clean))
    
    print("\n--- RESULTS ---")
    print(f"All {len(m_df)} matches:")
    print(f"  Median Delta V (DESI - SDSS): {med_all:+.3f} km/s")
    print(f"  Mean Delta V:                 {mean_all:+.3f} km/s (std: {std_all:.3f})")
    print(f"\nClean sample ({len(clean)} non-binary stars with |Delta V| < 15 km/s):")
    print(f"  Median Delta V (DESI - SDSS): {med_clean:+.3f} km/s")
    print(f"  Mean Delta V:                 {mean_clean:+.3f} +/- {sem_clean:.3f} km/s")
    print(f"  Standard Deviation:           {std_clean:.3f} km/s")
