import pandas as pd
import os

report_dir = "/gpfs/home/hce24xau/scratch/gen_kmers/data/african_sydney/kraken_reports/t2t_mcf_05_depletion_reports"
meta_file = "/gpfs/home/hce24xau/scratch/gen_kmers/data/african_sydney/geoethnic_master_batch1.csv"

# Get report samples
report_samples = {
    f.replace(".krakenuniq.report", "")
    for f in os.listdir(report_dir)
    if f.endswith(".report")
}

# Load metadata
df = pd.read_csv(meta_file)
meta_samples = set(df.iloc[:, 0])  # adjust column if needed

# Compare
print("In reports not in metadata:")
print(report_samples - meta_samples)

print("\nIn metadata not in reports:")
print(meta_samples - report_samples)
