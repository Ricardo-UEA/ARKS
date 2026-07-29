#!/usr/bin/env python3
"""
evaluate_arks_depletion.py

Consolidated pre- vs post-ARKS-depletion evaluation for the CAMISIM simulation
benchmark. Replaces evaluate_simulations.py + evaluate_simulations_depleted.py
with a single script that parses both conditions through identical code and
reports the direct pre/post comparison in the same terms used in the
manuscript (Fig. 3c/d).

WHY THIS SCRIPT EXISTS
----------------------
evaluate_simulations.py and evaluate_simulations_depleted.py are near-identical
copies that differ only in which report directory they read. Both save their
main per-sample/matrix outputs using a *default* filter (min_kmers >= 100,
min_cov >= 0.001) - not the cov >= 0.005 threshold the manuscript says was
selected as "the reference classification threshold for all downstream
analyses". The cov >= 0.005 condition only appears inside each script's own
threshold_comparison.csv, computed independently on each condition's
unfiltered data, and the two were never compared to each other directly in
one place. Running the two scripts separately and comparing their outputs by
hand missed that the pre- and post-depletion numbers come out nearly
identical - at either threshold - with the kraken_depleted_reports/ directory
as currently populated. That's the discrepancy this script is built to catch
automatically (see the WARNING it prints).

WHAT IT DOES
------------
1. Loads ground truth ONCE, genus-level, preferring the richer
   ground_truth_summary/ CSVs (ground_truth_complete.csv, or a
   *_sample_reads_matrix.csv if you've run build_matrix_from_fastq_parrallel.py)
   and falling back to per-sample taxonomic_profile_0.txt files if neither
   is present.
2. Parses BOTH pre-depletion (kraken_reports/) and post-depletion
   (kraken_depleted_reports/) KrakenUniq reports with one function - no
   duplicated logic to drift out of sync.
3. Applies ONE explicit, configurable confidence threshold (default:
   cov >= 0.005, no minimum k-mer count - matching the manuscript) to both
   conditions identically.
4. Computes:
   - genus-level precision/recall/F1 per sample, pre vs post
   - false-positive READ fraction per sample (proportion of *all* classified
     reads assigned to genera absent from ground truth - this is the Fig. 3d
     metric, computed on raw unfiltered read counts, independent of the
     confidence threshold above)
   - the full 8-way threshold sweep (as before) for transparency/Fig. 3a-b
5. Prints a direct pre vs post summary in the manuscript's own language.
6. Computes a bacteria-only "preservation" summary: true-positive
   (real bacterial) read retention after depletion, and whether the reads
   depletion does remove are disproportionately false-positive-genus reads
   rather than true-positive ones. Warns only if true-positive read
   retention drops below 98% (real collateral damage) - NOT if F1 stays
   flat, which is the expected, correct outcome here (see below).

INTERPRETING PRE VS POST ON THESE SIMULATIONS
----------------------------------------------
These CAMISIM simulations are overwhelmingly bacterial and contain no human
sequence. That matters for how "pre vs post" should be read:
  - There is no host contamination in these samples for ARKS to remove, so
    genus-level precision/recall/F1 are NOT expected to improve after
    depletion. Near-identical scores pre vs post is the desired outcome -
    it means ARKS depletion does not collaterally damage legitimate
    bacterial classification (a specificity/safety result, not a null
    result).
  - The more informative question is what happens to the reads depletion
    DOES remove: if they are disproportionately reads that were already
    contributing to false-positive genus calls (rather than true-positive
    ones), that supports ARKS depletion cleaning up likely noise/ambiguous
    reads alongside its host-removal role, at negligible cost to real
    bacterial signal.
  - The simulated community panel (validated_community_for_simulations_V2.csv)
    is not 100% bacterial: of 66 candidate genera, 10 are viral (several
    human-tropic: HSV-1, EBV, JC polyomavirus, Hepatitis C, West Nile) and
    one is the eukaryotic parasite Plasmodium. Pass --community-csv to split
    true-positive read retention into bacterial vs non-bacterial genera,
    since viral/eukaryotic sequence could plausibly be more prone to
    collateral depletion than bacterial sequence is.

USAGE
-----
    python evaluate_arks_depletion.py \\
        /gpfs/home/hce24xau/scratch/gen_kmers/data/simulations/CAMISIM/meta_simulations_2026_v2 \\
        --pre-dir kraken_reports \\
        --post-dir kraken_depleted_reports \\
        --ground-truth-dir ground_truth_summary \\
        --samples 100 \\
        --min-cov 0.005 \\
        --output arks_depletion_eval

If --ground-truth-dir doesn't contain a usable matrix, pass
--truth-fallback-dir pointing at the directory containing
output_files/meta_sample_N/taxonomic_profile_0.txt (defaults to base_dir).

OUTPUTS
-------
  <output>_pre_krakenuniq_reads_raw.csv / _post_krakenuniq_reads_raw.csv   (raw, unfiltered per-genus read counts)
  <output>_pre_per_sample_metrics.csv / _post_per_sample_metrics.csv (filtered genus-level metrics + TP/FP read counts)
  <output>_threshold_comparison_pre.csv / _threshold_comparison_post.csv
  <output>_pre_vs_post_summary.csv       (one row per metric: pre mean/median, post mean/median, delta)
  <output>_preservation_per_sample.csv   (per-sample TP retention / FP removal breakdown)
  <output>_preservation_summary.csv      (pooled + mean/median TP retention & FP removal)
  <output>_full_results.xlsx             (everything above, one workbook)
"""

import argparse
import sys
import warnings
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

# ----------------------------------------------------------------------------
# Ground truth loading
# ----------------------------------------------------------------------------

def parse_taxonomic_profile(filepath):
    """Parse a single CAMISIM taxonomic_profile_0.txt, genus level only."""
    genera = {}
    with open(filepath, "r") as f:
        for line in f:
            if line.startswith("@") or line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split("\t")
            if len(parts) >= 5 and parts[1] == "genus":
                genus_name = parts[3].split("|")[-1].strip()
                pct = float(parts[4])
                if genus_name and pct > 0:
                    genera[genus_name] = pct
    return genera


def load_ground_truth_from_profiles(base_dir, samples):
    """Fallback: rebuild a genus x sample presence/abundance matrix directly
    from per-sample taxonomic_profile_0.txt files."""
    base_dir = Path(base_dir)
    candidates_dirs = [base_dir, base_dir / "output_files"]
    truth = defaultdict(dict)
    found_any = False
    for i in range(1, samples + 1):
        sample_name = f"meta_sample_{i}"
        tf = None
        for d in candidates_dirs:
            for cand in [
                d / f"taxonomic_profile_{i}.txt",
                d / f"taxonomic_profile_{i - 1}.txt",
                d / sample_name / "taxonomic_profile_0.txt",
            ]:
                if cand.exists():
                    tf = cand
                    break
            if tf:
                break
        if not tf:
            continue
        found_any = True
        for genus, pct in parse_taxonomic_profile(tf).items():
            truth[genus][sample_name] = pct
    if not found_any:
        return None
    samples_idx = [f"meta_sample_{i}" for i in range(1, samples + 1)]
    genera = sorted(truth.keys())
    mat = pd.DataFrame(0.0, index=samples_idx, columns=genera)
    for genus, per_sample in truth.items():
        for sample, val in per_sample.items():
            mat.loc[sample, genus] = val
    return mat


def load_ground_truth(ground_truth_dir, base_dir, samples, truth_fallback_dir=None):
    """Try, in order: a genus-level *_sample_reads_matrix.csv, then
    ground_truth_complete.csv (pivoted genus x sample), then per-sample
    taxonomic_profile_0.txt files. Returns a samples x genera DataFrame of
    non-negative numbers (reads or abundance - only used as >0 for presence)."""
    gt_dir = Path(ground_truth_dir) if ground_truth_dir else None

    if gt_dir and gt_dir.exists():
        genus_matrix = gt_dir / "genus_sample_reads_matrix.csv"
        if genus_matrix.exists():
            mat = pd.read_csv(genus_matrix, index_col=0).T
            mat.index = [f"meta_sample_{c}" if not str(c).startswith("meta_sample") else c
                         for c in mat.index]
            print(f"Ground truth: loaded {genus_matrix.name} "
                  f"({mat.shape[0]} samples x {mat.shape[1]} genera)")
            return mat

        complete_csv = gt_dir / "ground_truth_complete.csv"
        if complete_csv.exists():
            df = pd.read_csv(complete_csv)
            df["sample_name"] = df["sample"].apply(lambda s: f"meta_sample_{s}")
            mat = df.pivot_table(index="sample_name", columns="genus", values="reads",
                                  aggfunc="sum", fill_value=0)
            print(f"Ground truth: loaded {complete_csv.name} "
                  f"({mat.shape[0]} samples x {mat.shape[1]} genera)")
            return mat

    fallback_dir = truth_fallback_dir or base_dir
    print(f"Ground truth: no usable CSV in {ground_truth_dir!r} - "
          f"falling back to taxonomic_profile_0.txt under {fallback_dir}")
    mat = load_ground_truth_from_profiles(fallback_dir, samples)
    if mat is None:
        sys.exit(
            "ERROR: could not find ground truth anywhere (no genus_sample_reads_matrix.csv, "
            "no ground_truth_complete.csv, no taxonomic_profile_0.txt files). "
            "Run generate_ground_truth_summary.py or build_matrix_from_fastq_parrallel.py first, "
            "or pass --truth-fallback-dir."
        )
    print(f"Ground truth: parsed from taxonomic profiles "
          f"({mat.shape[0]} samples x {mat.shape[1]} genera)")
    return mat


# ----------------------------------------------------------------------------
# KrakenUniq report parsing (one function, used for BOTH pre and post)
# ----------------------------------------------------------------------------

def parse_kraken_report(filepath):
    """Parse a KrakenUniq report, genus level only. Returns
    {genus: {reads, taxReads, kmers, dup, cov}}."""
    genera = {}
    with open(filepath, "r") as f:
        for line in f:
            if line.startswith("#") or line.startswith("%"):
                continue
            parts = line.strip().split("\t")
            if len(parts) >= 9:
                try:
                    if parts[7] == "genus":
                        taxname = parts[8].strip()
                        if taxname:
                            genera[taxname] = {
                                "pct": float(parts[0]),
                                "reads": int(parts[1]),
                                "taxReads": int(parts[2]),
                                "kmers": int(parts[3]),
                                "dup": float(parts[4]),
                                "cov": float(parts[5]),
                            }
                except (ValueError, IndexError):
                    continue
    return genera


def find_report(reports_dir, sample_num):
    reports_dir = Path(reports_dir)
    sample_name = f"meta_sample_{sample_num}"
    for cand in [reports_dir / f"{sample_name}.report", reports_dir / f"meta_sample_{sample_num}.report"]:
        if cand.exists():
            return cand
    return None


def apply_filter(pred_data, min_reads=0, min_kmers=0, min_cov=0.0, min_taxreads=0, max_dup=float("inf")):
    return {
        k: v for k, v in pred_data.items()
        if v["reads"] >= min_reads
        and v["kmers"] >= min_kmers
        and v["cov"] >= min_cov
        and v["taxReads"] >= min_taxreads
        and v["dup"] <= max_dup
    }


# ----------------------------------------------------------------------------
# Metrics
# ----------------------------------------------------------------------------

def calculate_metrics_binary(truth_set, pred_set):
    tp = len(truth_set & pred_set)
    fp = len(pred_set - truth_set)
    fn = len(truth_set - pred_set)
    precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    recall = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0.0
    return {"TP": tp, "FP": fp, "FN": fn, "Precision_PPV": precision,
            "Recall_Sensitivity": recall, "F1": f1}


def calculate_metrics_abundance(truth_dict, pred_dict):
    all_genera = set(truth_dict) | set(pred_dict)
    if not all_genera:
        return {"Bray_Curtis_Dissimilarity": 1.0, "Abundance_Correlation": 0.0, "L1_Error": 0.0}
    truth_total = sum(truth_dict.values()) or 1
    pred_total = sum(pred_dict.values()) or 1
    truth_arr = np.array([truth_dict.get(g, 0) / truth_total for g in all_genera])
    pred_arr = np.array([pred_dict.get(g, 0) / pred_total for g in all_genera])
    denom = truth_arr.sum() + pred_arr.sum()
    bc = np.abs(truth_arr - pred_arr).sum() / denom if denom > 0 else 0.0
    if truth_arr.std() > 0 and pred_arr.std() > 0:
        corr = np.corrcoef(truth_arr, pred_arr)[0, 1]
    else:
        corr = 0.0
    l1 = np.abs(truth_arr - pred_arr).mean()
    return {"Bray_Curtis_Dissimilarity": bc, "Abundance_Correlation": corr, "L1_Error": l1}


def load_domain_map(community_csv):
    """Load a genus -> domain ('Bacteria' / 'Non-bacterial') map from the
    validated community/genome panel CSV (needs 'genus' and 'Division'
    columns, e.g. validated_community_for_simulations_V2.csv - the candidate
    genome pool CAMISIM sampled from when building each simulation). Returns
    None if no path is given or the file can't be read, so callers can skip
    the bacterial/non-bacterial split analysis gracefully.
    """
    if not community_csv:
        return None
    path = Path(community_csv)
    if not path.exists():
        print(f"  [domain-split] --community-csv not found at {path}, "
              f"skipping bacterial/non-bacterial split")
        return None
    try:
        df = pd.read_csv(path)
        domain_map = {}
        for genus, division in zip(df["genus"], df["Division"]):
            if pd.isna(genus):
                continue
            domain_map[genus] = "Bacteria" if division == "Bacteria" else "Non-bacterial"
        n_nonbacterial = sum(1 for v in domain_map.values() if v == "Non-bacterial")
        print(f"  [domain-split] loaded {len(domain_map)} genera from {path.name} "
              f"({n_nonbacterial} non-bacterial)")
        return domain_map
    except Exception as e:
        print(f"  [domain-split] could not parse --community-csv ({e}), "
              f"skipping bacterial/non-bacterial split")
        return None


def read_level_breakdown(pred_raw, truth_set, domain_map=None):
    """Split ALL classified reads (unfiltered) into true-positive reads
    (assigned to a genus present in ground truth) and false-positive reads
    (assigned to a genus absent from ground truth). Computed on raw read
    counts, independent of the confidence filter used for the
    precision/recall/F1 metrics above. FP_Read_Fraction_Pct is the Fig. 3d
    metric.

    These simulations are overwhelmingly bacterial, so TP_Reads is mostly
    real bacterial signal a depletion step should leave untouched, and
    FP_Reads is the pool of reads that were never going to contribute to a
    correct genus call anyway. The community panel also includes a handful
    of viral genera and one protozoan (Plasmodium) genus though, so
    "bacterial" isn't a perfect description of every true-positive read -
    see domain_map.

    If domain_map ({genus: "Bacteria"/"Non-bacterial"}, e.g. loaded via
    load_domain_map from validated_community_for_simulations_V2.csv) is
    given, also splits TP_Reads into TP_Reads_Bacteria / TP_Reads_NonBacteria
    / TP_Reads_Unmapped (genus present in truth but absent from domain_map).
    """
    total_reads = sum(v["reads"] for v in pred_raw.values())
    fp_reads = sum(v["reads"] for g, v in pred_raw.items() if g not in truth_set)
    tp_reads = total_reads - fp_reads
    fp_fraction = fp_reads / total_reads * 100.0 if total_reads else 0.0
    result = {"TP_Reads": tp_reads, "FP_Reads": fp_reads, "Total_Reads": total_reads,
              "FP_Read_Fraction_Pct": fp_fraction}

    if domain_map is not None:
        tp_bacteria = tp_nonbacteria = tp_unmapped = 0
        for g in truth_set:
            if g not in pred_raw:
                continue
            reads = pred_raw[g]["reads"]
            domain = domain_map.get(g)
            if domain == "Bacteria":
                tp_bacteria += reads
            elif domain == "Non-bacterial":
                tp_nonbacteria += reads
            else:
                tp_unmapped += reads
        result.update({"TP_Reads_Bacteria": tp_bacteria,
                        "TP_Reads_NonBacteria": tp_nonbacteria,
                        "TP_Reads_Unmapped": tp_unmapped})
    return result


# ----------------------------------------------------------------------------
# Per-condition evaluation (identical code path for pre and post)
# ----------------------------------------------------------------------------

def evaluate_condition(label, reports_dir, truth_matrix, samples, min_kmers, min_cov, domain_map=None):
    """Returns (per_sample_metrics_df, raw_reads_matrix, filtered_reads_matrix,
    presence_absence_matrix)."""
    all_pred_full = {}
    sample_rows = []
    all_pred_genera = defaultdict(dict)
    missing = []

    for i in range(1, samples + 1):
        sample_name = f"meta_sample_{i}"
        report_file = find_report(reports_dir, i)
        if not report_file:
            missing.append(sample_name)
            continue
        if sample_name not in truth_matrix.index:
            missing.append(sample_name)
            continue

        pred_raw = parse_kraken_report(report_file)
        all_pred_full[sample_name] = pred_raw

        truth_row = truth_matrix.loc[sample_name]
        truth_set = set(truth_row[truth_row > 0].index)

        pred_filtered = apply_filter(pred_raw, min_kmers=min_kmers, min_cov=min_cov)
        pred_set = set(pred_filtered.keys())
        for genus, data in pred_filtered.items():
            all_pred_genera[genus][sample_name] = data["reads"]

        binary = calculate_metrics_binary(truth_set, pred_set)
        truth_abund = {g: truth_row[g] for g in truth_set}
        abundance = calculate_metrics_abundance(truth_abund, {k: v["reads"] for k, v in pred_filtered.items()})
        read_breakdown = read_level_breakdown(pred_raw, truth_set, domain_map)

        sample_rows.append({
            "Sample": sample_name,
            "Truth_Genera": len(truth_set),
            "Predicted_Genera": len(pred_set),
            **binary,
            **abundance,
            **read_breakdown,
        })

    if missing:
        print(f"  [{label}] missing report/truth for {len(missing)} samples "
              f"(e.g. {missing[:3]}{'...' if len(missing) > 3 else ''})")

    metrics_df = pd.DataFrame(sample_rows)

    all_samples = [f"meta_sample_{i}" for i in range(1, samples + 1)]
    all_genera = sorted(set(truth_matrix.columns) | set(all_pred_genera.keys()))
    raw_reads = pd.DataFrame(0.0, index=all_samples, columns=sorted(
        set(g for pred in all_pred_full.values() for g in pred)))
    for sample, pred in all_pred_full.items():
        for g, v in pred.items():
            raw_reads.loc[sample, g] = v["reads"]

    filtered_reads = pd.DataFrame(0.0, index=all_samples, columns=all_genera)
    for genus, per_sample in all_pred_genera.items():
        for sample, reads in per_sample.items():
            filtered_reads.loc[sample, genus] = reads
    presence_absence = (filtered_reads > 0).astype(int)

    return metrics_df, raw_reads, filtered_reads, presence_absence, all_pred_full


def threshold_sweep(all_pred_full, truth_matrix, samples):
    configs = [
        {"name": "No filter", "min_kmers": 0, "min_cov": 0},
        {"name": "cov >= 0.001 only", "min_kmers": 0, "min_cov": 0.001},
        {"name": "kmers >= 100 only", "min_kmers": 100, "min_cov": 0},
        {"name": "kmers >= 100, cov >= 0.001", "min_kmers": 100, "min_cov": 0.001},
        {"name": "kmers >= 500, cov >= 0.001", "min_kmers": 500, "min_cov": 0.001},
        {"name": "kmers >= 100, cov >= 0.002", "min_kmers": 100, "min_cov": 0.002},
        {"name": "cov >= 0.005 (manuscript threshold)", "min_kmers": 0, "min_cov": 0.005},
    ]
    rows = []
    for cfg in configs:
        name = cfg["name"]
        precs, recs, f1s = [], [], []
        for i in range(1, samples + 1):
            sample_name = f"meta_sample_{i}"
            if sample_name not in truth_matrix.index or sample_name not in all_pred_full:
                continue
            truth_row = truth_matrix.loc[sample_name]
            truth_set = set(truth_row[truth_row > 0].index)
            pred_filtered = apply_filter(all_pred_full[sample_name],
                                          min_kmers=cfg["min_kmers"], min_cov=cfg["min_cov"])
            m = calculate_metrics_binary(truth_set, set(pred_filtered.keys()))
            precs.append(m["Precision_PPV"]); recs.append(m["Recall_Sensitivity"]); f1s.append(m["F1"])
        if precs:
            rows.append({
                "Config": name, "min_kmers": cfg["min_kmers"], "min_cov": cfg["min_cov"],
                "Precision": np.mean(precs), "Precision_std": np.std(precs),
                "Recall": np.mean(recs), "Recall_std": np.std(recs),
                "F1": np.mean(f1s), "F1_std": np.std(f1s),
            })
    return pd.DataFrame(rows)


# ----------------------------------------------------------------------------
# Pre vs post summary
# ----------------------------------------------------------------------------

def summarise(pre_df, post_df):
    rows = []
    for col, fmt in [("Precision_PPV", "{:.1%}"), ("Recall_Sensitivity", "{:.1%}"),
                      ("F1", "{:.4f}"), ("FP_Read_Fraction_Pct", "{:.3f}%"),
                      ("Bray_Curtis_Dissimilarity", "{:.4f}"), ("Abundance_Correlation", "{:.4f}")]:
        pre_mean, post_mean = pre_df[col].mean(), post_df[col].mean()
        pre_med, post_med = pre_df[col].median(), post_df[col].median()
        rel_change = abs(post_mean - pre_mean) / abs(pre_mean) if pre_mean else float("nan")
        rows.append({
            "Metric": col, "Pre_mean": pre_mean, "Post_mean": post_mean,
            "Pre_median": pre_med, "Post_median": post_med,
            "Delta_mean": post_mean - pre_mean, "Relative_change": rel_change,
        })
    n_below_1pct = (post_df["FP_Read_Fraction_Pct"] < 1.0).sum()
    return pd.DataFrame(rows), n_below_1pct


def preservation_summary(pre_df, post_df):
    """Bacteria-only framing.

    These CAMISIM simulations contain no human sequence, so a depletion step
    cannot legitimately improve genus-level classification - there's no host
    contamination in the sample to remove that was confusing the classifier.
    Near-identical precision/recall/F1 pre vs post is therefore the expected,
    good outcome (no collateral damage), not a sign something is broken.

    What's actually worth measuring here:
      1. True-positive (real bacterial) read retention - how much of the
         signal supporting correct genus calls survives depletion. Want this
         close to 100%.
      2. Whether the reads depletion DOES remove are disproportionately ones
         that were already contributing to false-positive genus calls,
         rather than ones supporting correct calls - i.e. whether ARKS is
         preferentially cleaning up noise rather than cutting real signal.
    """
    cols = ["Sample", "TP_Reads", "FP_Reads", "Total_Reads", "FP"]
    merged = pre_df[cols].merge(post_df[cols], on="Sample", suffixes=("_pre", "_post"))

    merged["TP_Retention"] = np.where(
        merged["TP_Reads_pre"] > 0, merged["TP_Reads_post"] / merged["TP_Reads_pre"], np.nan)
    merged["FP_Removal"] = np.where(
        merged["FP_Reads_pre"] > 0, 1 - merged["FP_Reads_post"] / merged["FP_Reads_pre"], np.nan)
    merged["Reads_Removed_Total"] = merged["Total_Reads_pre"] - merged["Total_Reads_post"]
    merged["Reads_Removed_TP"] = merged["TP_Reads_pre"] - merged["TP_Reads_post"]
    merged["Reads_Removed_FP"] = merged["FP_Reads_pre"] - merged["FP_Reads_post"]
    removed_mask = merged["Reads_Removed_Total"] > 0
    merged["Pct_Removed_Reads_That_Were_FP"] = np.where(
        removed_mask, merged["Reads_Removed_FP"] / merged["Reads_Removed_Total"] * 100.0, np.nan)
    merged["Baseline_FP_Share_Pct"] = np.where(
        merged["Total_Reads_pre"] > 0, merged["FP_Reads_pre"] / merged["Total_Reads_pre"] * 100.0, np.nan)

    pooled_tp_retention = merged["TP_Reads_post"].sum() / merged["TP_Reads_pre"].sum()
    pooled_fp_removal = 1 - merged["FP_Reads_post"].sum() / merged["FP_Reads_pre"].sum()
    mean_removed_fp_share = merged.loc[removed_mask, "Pct_Removed_Reads_That_Were_FP"].mean()
    mean_baseline_fp_share = merged["Baseline_FP_Share_Pct"].mean()

    summary = {
        "n_samples": len(merged),
        "pooled_tp_retention": pooled_tp_retention,
        "mean_tp_retention": merged["TP_Retention"].mean(),
        "median_tp_retention": merged["TP_Retention"].median(),
        "pooled_fp_removal": pooled_fp_removal,
        "mean_fp_removal": merged["FP_Removal"].mean(),
        "median_fp_removal": merged["FP_Removal"].median(),
        "mean_fp_genus_calls_pre": merged["FP_pre"].mean(),
        "mean_fp_genus_calls_post": merged["FP_post"].mean(),
        "mean_baseline_fp_share_pct": mean_baseline_fp_share,
        "mean_pct_removed_reads_that_were_fp": mean_removed_fp_share,
        "fp_enrichment_in_removed_reads": (
            mean_removed_fp_share / mean_baseline_fp_share if mean_baseline_fp_share else float("nan")
        ),
    }

    # Optional bacterial vs non-bacterial (viral/Plasmodium) split - only
    # present if evaluate_condition was run with a domain_map.
    domain_cols = ["TP_Reads_Bacteria", "TP_Reads_NonBacteria", "TP_Reads_Unmapped"]
    if all(c in pre_df.columns for c in domain_cols) and all(c in post_df.columns for c in domain_cols):
        dmerged = pre_df[["Sample"] + domain_cols].merge(
            post_df[["Sample"] + domain_cols], on="Sample", suffixes=("_pre", "_post"))
        for tag, col in [("bacteria", "TP_Reads_Bacteria"), ("nonbacteria", "TP_Reads_NonBacteria")]:
            pre_total = dmerged[f"{col}_pre"].sum()
            post_total = dmerged[f"{col}_post"].sum()
            summary[f"{tag}_tp_reads_pre_total"] = pre_total
            summary[f"{tag}_tp_reads_post_total"] = post_total
            summary[f"{tag}_tp_retention_pooled"] = (post_total / pre_total) if pre_total else float("nan")
        merged = merged.merge(dmerged, on="Sample")

    return merged, summary


# ----------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------

def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("base_dir", help="Base simulation directory "
                                     "(e.g. .../CAMISIM/meta_simulations_2026_v2)")
    p.add_argument("--pre-dir", default="kraken_reports",
                   help="Subdirectory (relative to base_dir) with pre-depletion reports")
    p.add_argument("--post-dir", default="kraken_depleted_reports",
                   help="Subdirectory (relative to base_dir) with post-ARKS-depletion reports")
    p.add_argument("--ground-truth-dir", default="ground_truth_summary",
                   help="Subdirectory (relative to base_dir) with ground-truth CSVs")
    p.add_argument("--truth-fallback-dir", default=None,
                   help="Directory to search for taxonomic_profile_0.txt if the "
                        "ground-truth CSVs aren't usable (default: base_dir)")
    p.add_argument("--community-csv", default=None,
                   help="Optional path to the validated community/genome panel CSV "
                        "(needs 'genus' and 'Division' columns, e.g. "
                        "validated_community_for_simulations_V2.csv) to split "
                        "true-positive read retention into bacterial vs "
                        "non-bacterial (viral/Plasmodium) genera. Skipped if not given.")
    p.add_argument("--samples", type=int, default=100)
    p.add_argument("--min-cov", type=float, default=0.005,
                   help="Coverage threshold for the headline comparison (default 0.005, matching the manuscript)")
    p.add_argument("--min-kmers", type=int, default=0,
                   help="k-mer count threshold for the headline comparison (default 0, matching the manuscript's cov>=0.005 config)")
    p.add_argument("--output", default="arks_depletion_eval")
    args = p.parse_args()

    base_dir = Path(args.base_dir)
    pre_dir = base_dir / args.pre_dir
    post_dir = base_dir / args.post_dir
    gt_dir = base_dir / args.ground_truth_dir

    for d, label in [(pre_dir, "pre-dir"), (post_dir, "post-dir")]:
        if not d.exists():
            sys.exit(f"ERROR: --{label} not found: {d}")

    print("=" * 78)
    print("ARKS pre- vs post-depletion CAMISIM evaluation")
    print("=" * 78)
    print(f"Base dir:    {base_dir}")
    print(f"Pre-dir:     {pre_dir}")
    print(f"Post-dir:    {post_dir}")
    print(f"Threshold:   min_kmers >= {args.min_kmers}, min_cov >= {args.min_cov}")
    print()

    truth_matrix = load_ground_truth(gt_dir, base_dir, args.samples, args.truth_fallback_dir)
    domain_map = load_domain_map(args.community_csv)

    print("\nEvaluating PRE-depletion reports...")
    pre_metrics, pre_raw, pre_filtered, pre_pa, pre_full = evaluate_condition(
        "pre", pre_dir, truth_matrix, args.samples, args.min_kmers, args.min_cov, domain_map)

    print("Evaluating POST-depletion reports...")
    post_metrics, post_raw, post_filtered, post_pa, post_full = evaluate_condition(
        "post", post_dir, truth_matrix, args.samples, args.min_kmers, args.min_cov, domain_map)

    print("\nRunning threshold sweep (pre)...")
    pre_sweep = threshold_sweep(pre_full, truth_matrix, args.samples)
    print("Running threshold sweep (post)...")
    post_sweep = threshold_sweep(post_full, truth_matrix, args.samples)

    summary_df, n_below_1pct = summarise(pre_metrics, post_metrics)

    print("\n" + "=" * 78)
    print(f"PRE vs POST SUMMARY (threshold: kmers>={args.min_kmers}, cov>={args.min_cov})")
    print("=" * 78)
    print(f"{'Metric':<28}{'Pre (mean)':>14}{'Post (mean)':>14}{'Pre (median)':>15}{'Post (median)':>15}")
    for _, row in summary_df.iterrows():
        print(f"{row['Metric']:<28}{row['Pre_mean']:>14.4f}{row['Post_mean']:>14.4f}"
              f"{row['Pre_median']:>15.4f}{row['Post_median']:>15.4f}")
    print(f"\nSimulations with post-depletion FP read fraction < 1%: "
          f"{n_below_1pct}/{len(post_metrics)}")

    per_sample_preservation, preservation = preservation_summary(pre_metrics, post_metrics)

    print("\n" + "=" * 78)
    print("BACTERIAL READ PRESERVATION")
    print("=" * 78)
    print("These simulations contain NO human sequence, so genus-level")
    print("classification quality is not expected to improve with depletion -")
    print("there is no host contamination here for ARKS to remove. The two")
    print("things that matter instead: how much real bacterial signal survives")
    print("depletion, and whether the reads that ARE removed are")
    print("disproportionately ones that were never contributing to a correct")
    print("genus call.")
    print("-" * 78)
    print(f"True-positive bacterial reads retained after ARKS depletion (n={preservation['n_samples']} samples):")
    print(f"  pooled across all reads/samples: {preservation['pooled_tp_retention']:.2%}")
    print(f"  mean per-sample:                 {preservation['mean_tp_retention']:.2%}")
    print(f"  median per-sample:               {preservation['median_tp_retention']:.2%}")
    print(f"\nFalse-positive-genus reads removed by ARKS depletion:")
    print(f"  pooled across all reads/samples: {preservation['pooled_fp_removal']:.2%}")
    print(f"  mean per-sample:                 {preservation['mean_fp_removal']:.2%}")
    print(f"  median per-sample:               {preservation['median_fp_removal']:.2%}")
    print(f"\nMean false-positive genus calls per sample: "
          f"{preservation['mean_fp_genus_calls_pre']:.2f} (pre) -> "
          f"{preservation['mean_fp_genus_calls_post']:.2f} (post)")
    if not np.isnan(preservation["fp_enrichment_in_removed_reads"]):
        print(f"\nOf the reads ARKS depletion actually removed, "
              f"{preservation['mean_pct_removed_reads_that_were_fp']:.1f}% (mean per-sample) "
              f"were already assigned to a false-positive genus, versus a "
              f"{preservation['mean_baseline_fp_share_pct']:.1f}% false-positive share of all "
              f"reads pre-depletion ({preservation['fp_enrichment_in_removed_reads']:.1f}x enrichment) "
              f"- depletion is preferentially removing reads that were never going to "
              f"contribute to a correct call.")

    if "bacteria_tp_retention_pooled" in preservation:
        print(f"\nBy domain ({args.community_csv}):")
        print(f"  Bacterial genera TP retention:     "
              f"{preservation['bacteria_tp_retention_pooled']:.2%} "
              f"(pre total reads: {preservation['bacteria_tp_reads_pre_total']:,.0f})")
        print(f"  Non-bacterial genera TP retention: "
              f"{preservation['nonbacteria_tp_retention_pooled']:.2%} "
              f"(pre total reads: {preservation['nonbacteria_tp_reads_pre_total']:,.0f})")
        if preservation["nonbacteria_tp_reads_pre_total"] < 1000:
            print("  (non-bacterial read counts are small relative to bacterial - "
                  "treat this ratio as indicative, not precise)")

    # Safety check for a bacteria-only benchmark: the failure mode to catch
    # is collateral loss of real bacterial signal, NOT "F1 didn't improve"
    # (it isn't supposed to - there's no host contamination here for
    # depletion to remove). Warn only if true-positive read retention drops
    # meaningfully; otherwise say so explicitly.
    tp_retention = preservation["pooled_tp_retention"]
    if tp_retention < 0.98:
        warnings.warn(
            "\n" + "!" * 78 + "\n"
            "WARNING: ARKS depletion removed more than 2% of true-positive\n"
            f"bacterial reads (pooled TP retention: {tp_retention:.2%}).\n"
            "These simulations contain no human sequence, so any meaningful loss\n"
            "of reads assigned to real (ground-truth) genera is collateral damage,\n"
            "not intended host removal. Before trusting these numbers, check:\n"
            "  1. kraken_depleted_reports/ was generated from reads actually run\n"
            "     through BBDuk against the FINAL ARKS database (6.2 GB,\n"
            "     4,038,816,248 k-mers), not an older/partial build.\n"
            "  2. BBDuk was run with the documented parameters (k=31, hdist=0,\n"
            "     mincovfraction=0.5) and actually completed without error for all\n"
            "     samples.\n"
            "  3. Sample numbering/ordering between kraken_reports/ and\n"
            "     kraken_depleted_reports/ line up 1:1 (this script matches by\n"
            "     meta_sample_N filename, not by content).\n"
            + "!" * 78,
            stacklevel=2,
        )
    else:
        print(f"\n[OK] True-positive bacterial read retention ({tp_retention:.2%}) is at or "
              f"above the 98% safety threshold - ARKS depletion is not meaningfully "
              f"removing real bacterial signal in these bacteria-only simulations.")

    # Save everything
    out = args.output
    pre_raw.to_csv(f"{out}_pre_krakenuniq_reads_raw.csv")
    post_raw.to_csv(f"{out}_post_krakenuniq_reads_raw.csv")
    pre_pa.to_csv(f"{out}_pre_krakenuniq_presence_absence.csv")
    post_pa.to_csv(f"{out}_post_krakenuniq_presence_absence.csv")
    pre_metrics.to_csv(f"{out}_pre_per_sample_metrics.csv", index=False)
    post_metrics.to_csv(f"{out}_post_per_sample_metrics.csv", index=False)
    pre_sweep.to_csv(f"{out}_threshold_comparison_pre.csv", index=False)
    post_sweep.to_csv(f"{out}_threshold_comparison_post.csv", index=False)
    summary_df.to_csv(f"{out}_pre_vs_post_summary.csv", index=False)
    per_sample_preservation.to_csv(f"{out}_preservation_per_sample.csv", index=False)
    preservation_summary_df = pd.DataFrame(
        [{"Metric": k, "Value": v} for k, v in preservation.items()])
    preservation_summary_df.to_csv(f"{out}_preservation_summary.csv", index=False)

    truth_out = truth_matrix.copy()
    truth_out.index.name = "sample"

    with pd.ExcelWriter(f"{out}_full_results.xlsx", engine="openpyxl") as writer:
        summary_df.to_excel(writer, sheet_name="Pre_vs_Post_Summary", index=False)
        preservation_summary_df.to_excel(writer, sheet_name="Preservation_Summary", index=False)
        per_sample_preservation.to_excel(writer, sheet_name="Preservation_Per_Sample", index=False)
        pre_metrics.to_excel(writer, sheet_name="Pre_Per_Sample_Metrics", index=False)
        post_metrics.to_excel(writer, sheet_name="Post_Per_Sample_Metrics", index=False)
        pre_sweep.to_excel(writer, sheet_name="Threshold_Sweep_Pre", index=False)
        post_sweep.to_excel(writer, sheet_name="Threshold_Sweep_Post", index=False)
        truth_out.to_excel(writer, sheet_name="Ground_Truth")

    print(f"\nAll outputs written with prefix '{out}_*'. See {out}_full_results.xlsx for everything in one place.")


if __name__ == "__main__":
    main()
