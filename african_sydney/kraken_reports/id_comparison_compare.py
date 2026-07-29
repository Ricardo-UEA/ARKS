import pandas as pd

a = pd.read_csv("/gpfs/home/hce24xau/scratch/gen_kmers/data/african_sydney/kraken_reports/heroic1k_results/heroic1k_results_no_depletion/meta_data_analysis_included/per_sample_summary.csv")
b = pd.read_csv("/gpfs/home/hce24xau/scratch/gen_kmers/data/african_sydney/kraken_reports/heroic1k_results/heroic1k_results_depletion_panhuman_mcf_05/per_sample_summary.csv")

def sample_col(df):
    for c in df.columns:
        if "sample" in c.lower():
            return c
    return df.columns[0]

ca = sample_col(a)
cb = sample_col(b)

sa = set(a[ca].astype(str))
sb = set(b[cb].astype(str))

print("no_dep column:", ca)
print("dep column:", cb)
print("n no_dep:", len(sa))
print("n dep:", len(sb))
print("common:", len(sa & sb))

print("\nOnly in no_dep (first 10):")
print(sorted(list(sa - sb))[:10])

print("\nOnly in dep (first 10):")
print(sorted(list(sb - sa))[:10])


import pandas as pd

no_dep = "/gpfs/home/hce24xau/scratch/gen_kmers/data/african_sydney/kraken_reports/heroic1k_results/heroic1k_results_no_depletion/meta_data_analysis_included/per_sample_summary.csv"
dep    = "/gpfs/home/hce24xau/scratch/gen_kmers/data/african_sydney/kraken_reports/heroic1k_results/heroic1k_results_depletion_panhuman_mcf_05/per_sample_summary.csv"

a = pd.read_csv(no_dep)
b = pd.read_csv(dep)

a["sample_id"] = a["sample_id"].astype(str).str.strip()
b["sample_id"] = b["sample_id"].astype(str).str.strip()

sa = set(a["sample_id"])
sb = set(b["sample_id"])

print("no_dep n:", len(sa))
print("dep n:", len(sb))
print("common:", len(sa & sb))
print("\nfirst 10 no_dep:", sorted(sa)[:10])
print("first 10 dep:", sorted(sb)[:10])
print("\nonly in no_dep:", sorted(sa - sb)[:20])
print("\nonly in dep:", sorted(sb - sa)[:20])

