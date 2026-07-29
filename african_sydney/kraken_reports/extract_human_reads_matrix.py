#!/usr/bin/env python3
"""
ARKS – Human Read Fraction Matrix Extractor
============================================
Extracts Homo sapiens read counts from KrakenUniq reports across four
conditions and outputs a tidy CSV matrix ready for plotting in R.

KrakenUniq report column order (confirmed from real reports):
  0: %
  1: reads        ← use this for human_reads and total_reads
  2: taxReads
  3: kmers
  4: dup
  5: cov
  6: taxID
  7: rank
  8: taxName      ← name is col 8 (may be indented)

Total reads = unclassified (taxID 0) + root (taxID 1)
Human reads = reads on Homo sapiens line (taxID 9606)

Sample ID extracted from filename: KAL0054.krakenuniq.report → KAL0054
"""

import os
import re
import glob
import pandas as pd

#  PATHS

DIRS = {
    "No depletion": "/gpfs/home/hce24xau/scratch/gen_kmers/data/african_sydney/kraken_reports/no_depletion_kraken_reports",
    "GRCh38":       "/gpfs/home/hce24xau/scratch/gen_kmers/data/african_sydney/kraken_reports/GRCh38_mcf_05_depletion_reports",
    "T2T-CHM13":    "/gpfs/home/hce24xau/scratch/gen_kmers/data/african_sydney/kraken_reports/t2t_mcf_05_depletion_reports",
    "ARKS":         "/gpfs/home/hce24xau/scratch/gen_kmers/data/african_sydney/kraken_reports/ARKS_mcf_05_depletion_reports",
}

METADATA_PATH   = "/gpfs/home/hce24xau/scratch/gen_kmers/data/african_sydney/geoethnic_master_batch1.csv"
SAMPLE_ID_COL   = "SAMPLE ID"
ANCESTRY_COL    = "Ancestry (PS)"

OUT_DIR  = "/gpfs/home/hce24xau/scratch/gen_kmers/data/african_sydney/plots/"
OUT_FILE = os.path.join(OUT_DIR, "human_reads_matrix.csv")

CONDITION_ORDER = ["No depletion", "GRCh38", "T2T-CHM13", "ARKS"]

#  PARSE ONE KRAKENUNIQ REPORT

def parse_kraken_report(path):
    """
    Confirmed column layout (KrakenUniq v0.5.7):
      0:%  1:reads  2:taxReads  3:kmers  4:dup  5:cov  6:taxID  7:rank  8:taxName

    Total reads = unclassified reads (taxID 0) + root reads (taxID 1)
    Human reads = reads col for taxID 9606
    """
    unclassified_reads = None
    root_reads         = None
    human_reads        = None

    try:
        with open(path, "r") as f:
            for line in f:
                line = line.rstrip("\n")

                # Skip comment/header lines
                if line.startswith("#") or line.strip() == "":
                    continue
                # Skip the column header line
                if line.startswith("%"):
                    continue

                parts = line.split("\t")
                if len(parts) < 9:
                    continue

                try:
                    taxid = parts[6].strip()
                    name  = parts[8].strip()
                    reads = int(parts[1].strip())
                except (ValueError, IndexError):
                    continue

                if taxid == "0" or name == "unclassified":
                    unclassified_reads = reads

                elif taxid == "1" or name == "root":
                    root_reads = reads

                elif taxid == "9606" or name == "Homo sapiens":
                    human_reads = reads

    except Exception as e:
        print(f"  [WARN] Could not parse {os.path.basename(path)}: {e}")
        return None

    # Total = unclassified + root
    if unclassified_reads is not None and root_reads is not None:
        total_reads = unclassified_reads + root_reads
    elif root_reads is not None:
        total_reads = root_reads
    else:
        print(f"  [WARN] No root/unclassified line in {os.path.basename(path)}")
        return None

    if human_reads is None:
        print(f"  [WARN] No Homo sapiens line (taxID 9606) in {os.path.basename(path)}")
        return None

    if total_reads == 0:
        print(f"  [WARN] Total reads = 0 in {os.path.basename(path)}")
        return None

    return {
        "human_reads": human_reads,
        "total_reads": total_reads,
        "pct_human":   round(human_reads / total_reads * 100, 6),
    }


def extract_sample_id(filepath):
    """
    KAL0054.krakenuniq.report  →  KAL0054
    5060_report.txt            →  5060
    Strips everything from the first dot or underscore onward.
    """
    name = os.path.basename(filepath)
    # Take everything before the first dot or underscore
    m = re.match(r"^([A-Za-z0-9]+)", name)
    return m.group(1) if m else name


#  EXTRACT

os.makedirs(OUT_DIR, exist_ok=True)

print("=" * 60)
print("ARKS – Human Read Fraction Extractor")
print("=" * 60)

records = []

for condition in CONDITION_ORDER:
    directory = DIRS[condition]

    if not os.path.isdir(directory):
        print(f"\n[ERROR] Directory not found: {directory}")
        continue

    # Match .krakenuniq.report and other common extensions
    files = []
    for pattern in ["*.krakenuniq.report", "*.report", "*.txt",
                    "*.tsv", "*.kreport", "*.kraken"]:
        files.extend(glob.glob(os.path.join(directory, pattern)))
    # Deduplicate and sort
    files = sorted(set(files))

    print(f"\n{condition}: {len(files)} files")
    if files:
        print(f"  Example: {os.path.basename(files[0])}")

    n_ok = 0
    for filepath in files:
        result = parse_kraken_report(filepath)
        if result is None:
            continue
        sample_id = extract_sample_id(filepath)
        records.append({
            "sample_id": sample_id,
            "condition": condition,
            **result,
        })
        n_ok += 1

    print(f"  Parsed OK: {n_ok}/{len(files)}")

if not records:
    raise RuntimeError(
        "\nNo data extracted. Check:\n"
        "  1. Directory paths exist on HPC\n"
        "  2. Files match *.krakenuniq.report\n"
        "  3. Run: head -20 <one_report> to inspect format\n"
    )

df = pd.DataFrame(records)
df["condition"] = pd.Categorical(df["condition"],
                                  categories=CONDITION_ORDER, ordered=True)

print(f"\n{'=' * 60}")
print(f"Extracted {len(df)} records | {df['sample_id'].nunique()} unique samples")

# ── Keep only samples present in ALL four conditions ──────────
counts       = df.groupby("sample_id")["condition"].nunique()
complete_ids = counts[counts == len(DIRS)].index.tolist()
incomplete   = counts[counts < len(DIRS)].index.tolist()

if incomplete:
    print(f"\n[INFO] {len(incomplete)} samples missing from ≥1 condition:")
    for s in incomplete[:10]:
        present = df[df["sample_id"] == s]["condition"].tolist()
        print(f"  {s}: present in {present}")
    if len(incomplete) > 10:
        print(f"  ... and {len(incomplete) - 10} more")

df = df[df["sample_id"].isin(complete_ids)].copy()
print(f"\n{len(complete_ids)} samples complete across all conditions")

# ── Join metadata ─────────────────────────────────────────────
if os.path.exists(METADATA_PATH):
    print(f"\nLoading metadata: {METADATA_PATH}")
    meta = pd.read_csv(METADATA_PATH)
    meta.columns = meta.columns.str.strip()

    sid  = SAMPLE_ID_COL.strip()
    anc  = ANCESTRY_COL.strip()

    if sid not in meta.columns:
        print(f"[WARN] '{sid}' not in metadata. Available: {list(meta.columns)[:8]}")
    elif anc not in meta.columns:
        print(f"[WARN] '{anc}' not in metadata. Available: {list(meta.columns)[:8]}")
    else:
        meta[sid] = meta[sid].astype(str).str.strip()
        df["sample_id"] = df["sample_id"].astype(str).str.strip()

        df = df.merge(
            meta[[sid, anc]].drop_duplicates(),
            left_on="sample_id", right_on=sid,
            how="left"
        )
        if sid != "sample_id":
            df = df.drop(columns=[sid])
        df = df.rename(columns={anc: "Ancestry"})

        matched   = df["Ancestry"].notna().sum()
        unmatched = df["Ancestry"].isna().sum()
        print(f"  Matched: {matched} | Unmatched: {unmatched}")
        if unmatched > 0:
            bad = df[df["Ancestry"].isna()]["sample_id"].unique()[:5]
            print(f"  Unmatched sample IDs (first 5): {bad}")
else:
    print(f"[WARN] Metadata not found — saving without ancestry")
    df["Ancestry"] = "Unknown"

# ── Summary ───────────────────────────────────────────────────
print(f"\n{'─' * 60}")
print("Median human read fraction (%) per condition:")
for cond in CONDITION_ORDER:
    sub = df[df["condition"] == cond]["pct_human"]
    if len(sub):
        print(f"  {cond:<16} n={len(sub):>4}  "
              f"median={sub.median():.4f}%  "
              f"IQR=[{sub.quantile(0.25):.4f}%, "
              f"{sub.quantile(0.75):.4f}%]")

# ── Save ──────────────────────────────────────────────────────
cols = ["sample_id", "condition", "human_reads",
        "total_reads", "pct_human", "Ancestry"]
cols = [c for c in cols if c in df.columns]
df[cols].to_csv(OUT_FILE, index=False)

print(f"\n{'=' * 60}")
print(f"✅ Saved: {OUT_FILE}")
print(f"   Rows: {len(df)} | Samples: {df['sample_id'].nunique()}")
print(f"\nNext: Rscript plot_human_reads.R")
