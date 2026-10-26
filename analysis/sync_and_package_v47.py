import os
import shutil
import zipfile

base_dir = "D:/ARGON"
artifact_dir = "C:/Users/PC/.gemini/antigravity/brain/10c6cb6d-2b0f-4971-9598-378151ec0b61"

arxiv_files = [
    "gd1_stream_kinematic_audit.tex",
    "gd1_stream_kinematic_audit.pdf",
    "v45_physical_attribution.png",
    "v46_orbital_phase_confirmation.png",
    "v47_kinematic_confrontation.png",
    "nibauer_bin_comparison.csv",
    "sdss_desi_empirical_crossmatches.csv"
]

zenodo_files = arxiv_files + [
    "zenodo.json",
    "measure_exact_zero_point.py",
    "compute_nibauer_bins_confrontation.py"
]

# Create arXiv zip
arxiv_zip = os.path.join(base_dir, "GD1_v47_arXiv_submission.zip")
with zipfile.ZipFile(arxiv_zip, "w", zipfile.ZIP_DEFLATED) as z:
    for f in arxiv_files:
        src = os.path.join(base_dir, f)
        z.write(src, arcname=f)
print(f"Created {arxiv_zip}")

# Create Zenodo zip
zenodo_zip = os.path.join(base_dir, "GD1_v47_Zenodo_deposit.zip")
with zipfile.ZipFile(zenodo_zip, "w", zipfile.ZIP_DEFLATED) as z:
    for f in zenodo_files:
        src = os.path.join(base_dir, f)
        z.write(src, arcname=f)
print(f"Created {zenodo_zip}")

# Copy to artifact directory
for f in ["gd1_stream_kinematic_audit.pdf", "GD1_v47_arXiv_submission.zip", "GD1_v47_Zenodo_deposit.zip", "v47_kinematic_confrontation.png", "v47_fundamental_physics_constraints.png"]:
    src = os.path.join(base_dir, f)
    dst = os.path.join(artifact_dir, f)
    shutil.copy2(src, dst)
    print(f"Copied {f} to artifacts directory.")
