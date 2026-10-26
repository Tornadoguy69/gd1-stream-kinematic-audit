import os
import urllib.request
import urllib.parse
import json
import numpy as np
from astropy.table import Table

print("=" * 80)
print("RESEARCH v46: MULTI-EPOCH RADIAL VELOCITY CROSS-MATCH")
print("=" * 80)

# Load Table7_all and v45 attribution summary
t = Table.read('D:/ARGON/Table7_all.fits').to_pandas()
t['SOURCE_ID'] = t['SOURCE_ID'].astype(str).str.strip()

with open('D:/ARGON/v45_results/v45_attribution_summary.json', 'r') as fp:
    v45 = json.load(fp)

tail_stars = v45['top_contributors']

print(f"Checking all {len(tail_stars)} continuous tail stars across SDSS DR17...")
sdss_matches = []

for s in tail_stars:
    sid = s['id']
    row = t[t['SOURCE_ID'] == sid].iloc[0]
    ra = float(row['RA'])
    dec = float(row['Dec'])
    v_desi = float(row['V_LOS'])
    ev_desi = float(row['V_ERR'])
    
    # 2.5 arcsecond cone search
    cosd = np.cos(np.radians(dec))
    r_deg = 2.5 / 3600.0
    
    sql = f"""
    SELECT s.specObjID, s.bestObjID, s.ra, s.dec, s.mjd, s.plate, s.fiberID,
           s.zwarning, s.z, s.zErr, s.class, s.subClass,
           (s.z * 299792.458) as v_rad_kms, (s.zErr * 299792.458) as e_vrad
    FROM specObjAll as s
    WHERE s.ra BETWEEN {ra - r_deg/cosd:.6f} AND {ra + r_deg/cosd:.6f}
      AND s.dec BETWEEN {dec - r_deg:.6f} AND {dec + r_deg:.6f}
      AND s.zwarning = 0
    """
    url = f"https://skyserver.sdss.org/dr17/SkyServerWS/SearchTools/SqlSearch?cmd={urllib.parse.quote(sql)}&format=json"
    try:
        req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64)'})
        with urllib.request.urlopen(req, timeout=15) as resp:
            data = json.loads(resp.read().decode('utf-8'))
            rows = data[0].get('Rows', [])
            if rows:
                for r in rows:
                    v_sdss = float(r['v_rad_kms'])
                    ev_sdss = float(r['e_vrad'])
                    mjd_sdss = float(r['mjd'])
                    plate = int(r['plate'])
                    fiber = int(r['fiberID'])
                    spec_id = str(r['specObjID'])
                    
                    delta_v = v_desi - v_sdss
                    sigma_delta_v = np.hypot(ev_desi, ev_sdss)
                    z_shift = abs(delta_v) / sigma_delta_v
                    
                    match_record = {
                        "rank": s['rank'],
                        "source_id": sid,
                        "classification": s['classification'],
                        "ra": ra,
                        "dec": dec,
                        "v_desi": v_desi,
                        "ev_desi": ev_desi,
                        "v_sdss": v_sdss,
                        "ev_sdss": ev_sdss,
                        "mjd_sdss": mjd_sdss,
                        "plate": plate,
                        "fiber": fiber,
                        "specObjID": spec_id,
                        "delta_v_kms": delta_v,
                        "sigma_delta_v": sigma_delta_v,
                        "z_shift": z_shift
                    }
                    sdss_matches.append(match_record)
                    print(f"MATCH: Rank {s['rank']:2d} ({sid}) | {s['classification']}")
                    print(f"   DESI  V_LOS = {v_desi:+7.2f} +/- {ev_desi:4.2f} km/s")
                    print(f"   SDSS  V_LOS = {v_sdss:+7.2f} +/- {ev_sdss:4.2f} km/s (MJD {mjd_sdss:.0f}, Plate {plate})")
                    print(f"   -> Velocity Shift Delta V = {delta_v:+7.2f} +/- {sigma_delta_v:4.2f} km/s ({z_shift:.2f} sigma)\n")
    except Exception as e:
        print(f"Error querying {sid}: {e}")

print(f"\nTotal SDSS DR17 matches found: {len(sdss_matches)}")

# Save to JSON
os.makedirs("D:/ARGON/v46_results", exist_ok=True)
with open("D:/ARGON/v46_results/v46_multi_epoch_matches.json", "w") as fp:
    json.dump(sdss_matches, fp, indent=2)
print("Saved matches to D:/ARGON/v46_results/v46_multi_epoch_matches.json")
