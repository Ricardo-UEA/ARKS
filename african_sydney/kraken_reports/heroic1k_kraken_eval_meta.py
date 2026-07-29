#!/usr/bin/env python3
"""
HEROIC1K KrakenUniq Report Evaluation Pipeline — with Metadata
================================================================
Parses KrakenUniq reports and integrates GEOETHNIC_MASTER metadata
for ancestry-aware, clinically-contextualised analysis.

Usage:

#Non depletion usage
python heroic1k_kraken_eval_meta.py \
  --report_dir /gpfs/home/hce24xau/scratch/gen_kmers/data/african_sydney/kraken_reports/no_depletion_reports \
  --metadata /gpfs/home/hce24xau/scratch/gen_kmers/data/african_sydney/geoethnic_master_batch1.csv \
  --output_dir /gpfs/home/hce24xau/scratch/gen_kmers/data/african_sydney/kraken_reports/heroic1k_results/heroic1k_results_no_depletion/meta_data_analysis_included \
  --cov_threshold 0.005 \
  --top_genera 40

#Panhuman depletion usage
python heroic1k_kraken_eval_meta.py \
  --report_dir /gpfs/home/hce24xau/scratch/gen_kmers/data/african_sydney/kraken_reports/panhuman_depletion_reports \
  --metadata /gpfs/home/hce24xau/scratch/gen_kmers/data/african_sydney/geoethnic_master_batch1.csv \
  --output_dir /gpfs/home/hce24xau/scratch/gen_kmers/data/african_sydney/kraken_reports/heroic1k_results/heroic1k_results_depletion_panhuman \
  --cov_threshold 0.005 \
  --top_genera 40 

#T2T depletion Evaluation 
python heroic1k_kraken_eval_meta.py \
  --report_dir /gpfs/home/hce24xau/scratch/gen_kmers/data/african_sydney/kraken_reports/t2t_depletion_reports \
  --metadata /gpfs/home/hce24xau/scratch/gen_kmers/data/african_sydney/geoethnic_master_batch1.csv \
  --output_dir /gpfs/home/hce24xau/scratch/gen_kmers/data/african_sydney/kraken_reports/heroic1k_results/heroic1k_results_depletion_t2t \
  --cov_threshold 0.005 \
  --top_genera 40 
"""

import argparse
import sys
from pathlib import Path
from collections import defaultdict

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
from scipy.stats import kruskal, spearmanr
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import pdist

def safe_kruskal(*groups):
    """Kruskal-Wallis only if >=2 groups have >1 sample and >1 unique value."""
    valid = [g for g in groups if len(g) > 1 and len(np.unique(g)) > 1]
    if len(valid) < 2:
        return None, None
    try:
        return kruskal(*valid)
    except ValueError:
        return None, None

# ── Constants ─────────────────────────────────────────────────────────────────

HUMAN_NAMES  = {"Homo sapiens", "Homo"}
HUMAN_TAXIDS = {"9606", "9605"}
CONTAMINANTS = {"Buchnera", "Lymphocryptovirus", "Epstein-Barr virus"}

# Ancestry colour palette — consistent across all plots
ANCESTRY_COLORS = {
    "Southern African": "#E53935",   # red
    "European":         "#1565C0",   # blue
    "Admixed":          "#6A1B9A",   # purple
    "Asian":            "#2E7D32",   # green
    "Unknown":          "#757575",   # grey
}

# ── Parser ─────────────────────────────────────────────────────────────────────

def parse_report(filepath):
    genera = []
    unclass_reads = 0
    root_reads    = 0
    human_reads   = 0

    with open(filepath) as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[0].strip() == "%":
                continue
            try:
                reads_clade  = int(parts[1].strip())
                reads_direct = int(parts[2].strip())
                kmers        = int(parts[3].strip())
                dup          = float(parts[4].strip())
                cov          = float(parts[5].strip())
                taxid        = parts[6].strip()
                rank         = parts[7].strip()
                name         = parts[8].strip()
            except (ValueError, IndexError):
                continue

            if name == "unclassified":
                unclass_reads = reads_clade
            if rank == "no rank" and name == "root":
                root_reads = reads_clade
            if name in HUMAN_NAMES or taxid in HUMAN_TAXIDS:
                if reads_clade > human_reads:
                    human_reads = reads_clade
            if rank == "genus":
                genera.append({
                    "reads_clade": reads_clade, "reads_direct": reads_direct,
                    "kmers": kmers, "dup": dup, "cov": cov,
                    "taxid": taxid, "name": name,
                })

    return genera, unclass_reads + root_reads, human_reads, unclass_reads, root_reads


def summarise_sample(sid, genera, total_reads, human_reads, unclass_reads, classified):
    total = max(total_reads, 1)
    return {
        "sample_id":           sid,
        "total_reads":         total_reads,
        "classified_reads":    classified,
        "unclassified_reads":  unclass_reads,
        "human_reads":         human_reads,
        "nonhuman_classified": classified - human_reads,
        "human_pct":           human_reads / total * 100,
        "unclassified_pct":    unclass_reads / total * 100,
        "classified_pct":      classified / total * 100,
        "n_genera_detected":   len(genera),
    }


# ── Diversity ─────────────────────────────────────────────────────────────────

def shannon(counts):
    a = np.array(counts, dtype=float); a = a[a > 0]
    if a.sum() == 0: return 0.0
    p = a / a.sum()
    return float(-np.sum(p * np.log(p)))

def simpson(counts):
    a = np.array(counts, dtype=float); n = a.sum()
    if n < 2: return 0.0
    return float(1 - np.sum(a * (a - 1)) / (n * (n - 1)))


# ── Metadata loader ───────────────────────────────────────────────────────────

def load_metadata(metadata_path):
    df = pd.read_csv(metadata_path, dtype=str)
    df.columns = df.columns.str.strip()
    df["SAMPLE ID"] = df["SAMPLE ID"].astype(str).str.strip()

    # Normalise key columns
    df["ancestry"]   = df["Ancestry (PS)"].str.strip().fillna("Unknown")
    df["broad_anc"]  = df["Broad_Ancestry (PS)"].str.strip().fillna("Unknown")
    df["ethnicity"]  = df["Ethnicity_Self_Rep"].str.strip().fillna("Unknown")
    df["country"]    = df["Country Recruitment (born)"].str.strip().fillna("Unknown")
    df["region"]     = df["Region/City Recruited"].str.strip().fillna("Unknown")
    df["isup"]       = pd.to_numeric(df["ISUP FINAL (VH)"], errors="coerce")
    df["age"]        = pd.to_numeric(df["Age_Diag"],         errors="coerce")
    df["psa"]        = pd.to_numeric(df["PSA"],              errors="coerce")

    return df.set_index("SAMPLE ID")


# ── Load all reports ───────────────────────────────────────────────────────────

def load_all(report_dir, pattern, cov_threshold):
    report_dir = Path(report_dir)
    files = sorted(report_dir.glob(pattern)) or sorted(report_dir.glob("*.report"))
    if not files:
        sys.exit(f"[ERROR] No report files found in {report_dir}")
    print(f"[INFO] Found {len(files)} report files")

    summaries   = []
    diversity   = []
    flagged     = []
    kmer_matrix = defaultdict(dict)

    for f in files:
        sid = f.name.split(".krakenuniq")[0].split(".report")[0]
        genera, total, human, unclass, classified = parse_report(f)
        summaries.append(summarise_sample(sid, genera, total, human, unclass, classified))

        kmer_counts = [g["kmers"] for g in genera if g["kmers"] > 0]
        diversity.append({
            "sample_id": sid, "observed_genera": len(genera),
            "shannon": shannon(kmer_counts), "simpson": simpson(kmer_counts),
        })

        for g in genera:
            kmer_matrix[g["name"]][sid] = g["kmers"]
            if g["cov"] >= cov_threshold:
                flagged.append({
                    "sample_id": sid, "genus": g["name"],
                    "reads_clade": g["reads_clade"], "kmers": g["kmers"],
                    "dup": g["dup"], "cov": g["cov"],
                    "is_contaminant": any(c.lower() in g["name"].lower() for c in CONTAMINANTS),
                })

    df_summ = pd.DataFrame(summaries)
    df_div  = pd.DataFrame(diversity)
    df_flag = pd.DataFrame(flagged) if flagged else pd.DataFrame(
        columns=["sample_id","genus","reads_clade","kmers","dup","cov","is_contaminant"])

    all_samples = [s["sample_id"] for s in summaries]
    df_kmer = pd.DataFrame(kmer_matrix, index=all_samples).fillna(0).astype(int)

    return df_summ, df_div, df_flag, df_kmer


# ── Merge with metadata ────────────────────────────────────────────────────────

def merge_metadata(df, meta):
    """Left-join sample table with metadata on sample_id."""
    merged = df.set_index("sample_id").join(meta, how="left")
    merged["ancestry"] = merged["ancestry"].fillna("Unknown")
    merged["country"]  = merged["country"].fillna("Unknown")
    n_matched = merged["ancestry"].ne("Unknown").sum()
    print(f"[INFO] Metadata matched: {n_matched}/{len(merged)} samples")
    return merged.reset_index().rename(columns={"index": "sample_id"})


# ── Helper: ancestry colour list ──────────────────────────────────────────────

def anc_colors(series):
    return [ANCESTRY_COLORS.get(a, "#757575") for a in series]


# ── Plots ──────────────────────────────────────────────────────────────────────

def plot_read_overview(df, out_dir):
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    fig.suptitle("HEROIC1K — Read-Level Overview (No-Depletion Baseline)",
                 fontsize=13, fontweight="bold")
    configs = [
        ("human_pct",        "Human Reads (%)",        "#E53935"),
        ("unclassified_pct", "Unclassified Reads (%)", "#9E9E9E"),
        ("total_reads",      "Total Reads (M)",         "#0288D1"),
    ]
    for ax, (col, label, color) in zip(axes, configs):
        vals = df[col] / (1e6 if col == "total_reads" else 1)
        ax.hist(vals, bins=30, color=color, edgecolor="white", linewidth=0.4)
        med = vals.median()
        ax.axvline(med, color="black", ls="--", lw=1.2, label=f"Median: {med:.2f}")
        ax.set_xlabel(label); ax.set_ylabel("Samples"); ax.set_title(label)
        ax.legend(fontsize=9)
    plt.tight_layout()
    plt.savefig(out_dir / "01_read_overview.png", dpi=150, bbox_inches="tight")
    plt.close(); print("[PLOT] 01_read_overview.png")


def plot_composition_bar(df, out_dir):
    df_s  = df.sort_values("human_pct", ascending=False).reset_index(drop=True)
    total = df_s["total_reads"].replace(0, 1)
    h  = df_s["human_reads"] / total * 100
    nc = df_s["nonhuman_classified"] / total * 100
    u  = df_s["unclassified_reads"] / total * 100

    fig, ax = plt.subplots(figsize=(max(14, len(df_s) * 0.1), 5))
    x = np.arange(len(df_s))
    ax.bar(x, h,            color="#E53935", width=1.0, label="Human")
    ax.bar(x, nc, bottom=h, color="#43A047", width=1.0, label="Non-human classified")
    ax.bar(x, u,  bottom=h+nc, color="#9E9E9E", width=1.0, label="Unclassified")
    ax.set_xlim(-0.5, len(df_s)-0.5); ax.set_ylim(0, 100)
    ax.set_xlabel(f"Samples (n={len(df_s)}, sorted by human %)"); ax.set_ylabel("Read Composition (%)")
    ax.set_title("HEROIC1K — Read Composition per Sample"); ax.set_xticks([]); ax.legend(loc="upper right")
    plt.tight_layout()
    plt.savefig(out_dir / "02_read_composition.png", dpi=150, bbox_inches="tight")
    plt.close(); print("[PLOT] 02_read_composition.png")


def plot_ancestry_read_stats(df, out_dir):
    """Boxplots of human%, unclassified%, total reads — stratified by ancestry."""
    anc_order = ["Southern African", "European", "Admixed", "Asian", "Unknown"]
    anc_order = [a for a in anc_order if a in df["ancestry"].values]

    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    fig.suptitle("HEROIC1K — Read-Level Metrics by Ancestry", fontsize=13, fontweight="bold")

    for ax, col, label in zip(
        axes,
        ["human_pct", "unclassified_pct", "total_reads"],
        ["Human Reads (%)", "Unclassified Reads (%)", "Total Reads"],
    ):
        data_by_anc = [df.loc[df["ancestry"] == a, col].dropna().values for a in anc_order]
        bp = ax.boxplot(data_by_anc, patch_artist=True, notch=False,
                        medianprops=dict(color="black", lw=2))
        for patch, anc in zip(bp["boxes"], anc_order):
            patch.set_facecolor(ANCESTRY_COLORS.get(anc, "#757575"))
            patch.set_alpha(0.8)
        ax.set_xticklabels(anc_order, rotation=20, ha="right", fontsize=9)
        ax.set_ylabel(label); ax.set_title(label)

        # Kruskal-Wallis p-value
        valid = [d for d in data_by_anc if len(d) > 1]
        if len(valid) >= 2:
            stat, p = safe_kruskal(*valid)
            if stat is None: continue
            ax.set_title(f"{label}\nKruskal-Wallis p={p:.3f}" if p is not None else label, fontsize=9)

    plt.tight_layout()
    plt.savefig(out_dir / "03_ancestry_read_stats.png", dpi=150, bbox_inches="tight")
    plt.close(); print("[PLOT] 03_ancestry_read_stats.png")


def plot_alpha_diversity(df, out_dir):
    """Alpha diversity histograms coloured by ancestry."""
    fig, axes = plt.subplots(1, 3, figsize=(15, 4))
    fig.suptitle("HEROIC1K — Genus-Level Alpha Diversity", fontsize=12, fontweight="bold")

    anc_order = ["Southern African", "European", "Admixed", "Asian", "Unknown"]
    anc_present = [a for a in anc_order if a in df["ancestry"].values]

    for ax, col, label, color in zip(
        axes,
        ["observed_genera", "shannon", "simpson"],
        ["Observed Genera", "Shannon Entropy (k-mers)", "Simpson Index (k-mers)"],
        ["#7B1FA2", "#0288D1", "#F4511E"],
    ):
        ax.hist(df[col], bins=25, color=color, edgecolor="white", linewidth=0.4)
        med = df[col].median()
        ax.axvline(med, color="black", ls="--", lw=1.2, label=f"Median: {med:.2f}")
        ax.set_xlabel(label); ax.set_ylabel("Samples"); ax.set_title(label)
        ax.legend(fontsize=9)
    plt.tight_layout()
    plt.savefig(out_dir / "04_alpha_diversity.png", dpi=150, bbox_inches="tight")
    plt.close(); print("[PLOT] 04_alpha_diversity.png")


def plot_ancestry_diversity(df, out_dir):
    """Shannon and observed genera by ancestry — boxplots + strip."""
    anc_order = ["Southern African", "European", "Admixed", "Asian", "Unknown"]
    anc_order = [a for a in anc_order if a in df["ancestry"].values]

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    fig.suptitle("HEROIC1K — Alpha Diversity by Ancestry", fontsize=12, fontweight="bold")

    for ax, col, label in zip(
        axes,
        ["shannon", "observed_genera"],
        ["Shannon Entropy", "Observed Genera"],
    ):
        data_by_anc = [df.loc[df["ancestry"] == a, col].dropna().values for a in anc_order]
        bp = ax.boxplot(data_by_anc, patch_artist=True,
                        medianprops=dict(color="black", lw=2), zorder=2)
        for patch, anc in zip(bp["boxes"], anc_order):
            patch.set_facecolor(ANCESTRY_COLORS.get(anc, "#757575")); patch.set_alpha(0.7)

        # Overlay strip
        for i, (anc, vals) in enumerate(zip(anc_order, data_by_anc)):
            jitter = np.random.default_rng(42).uniform(-0.15, 0.15, len(vals))
            ax.scatter(i + 1 + jitter, vals,
                       color=ANCESTRY_COLORS.get(anc, "#757575"),
                       alpha=0.5, s=15, zorder=3)

        valid = [d for d in data_by_anc if len(d) > 1]
        if len(valid) >= 2:
            stat, p = safe_kruskal(*valid)
            if stat is None: continue
            ax.set_title(f"{label}\nKruskal-Wallis p={p:.3f}" if p is not None else label, fontsize=9)
        else:
            ax.set_title(label)

        ax.set_xticklabels(anc_order, rotation=20, ha="right", fontsize=9)
        ax.set_ylabel(label)

    plt.tight_layout()
    plt.savefig(out_dir / "05_ancestry_diversity.png", dpi=150, bbox_inches="tight")
    plt.close(); print("[PLOT] 05_ancestry_diversity.png")


def plot_depth_diversity_scatter(df, out_dir):
    """Scatter: total reads vs Shannon — coloured by ancestry, with Spearman r."""
    fig, ax = plt.subplots(figsize=(8, 6))

    anc_order = ["Southern African", "European", "Admixed", "Asian", "Unknown"]
    for anc in anc_order:
        sub = df[df["ancestry"] == anc]
        if sub.empty: continue
        ax.scatter(sub["total_reads"] / 1e6, sub["shannon"],
                   color=ANCESTRY_COLORS.get(anc, "#757575"),
                   label=f"{anc} (n={len(sub)})", alpha=0.7, s=35, edgecolors="white", lw=0.3)

    r, p = spearmanr(df["total_reads"], df["shannon"])
    ax.set_xlabel("Total Reads (M)"); ax.set_ylabel("Shannon Entropy")
    ax.set_title(f"Read Depth vs Diversity\nSpearman r={r:.2f}, p={p:.3e}", fontsize=11)
    ax.legend(fontsize=8, loc="lower right")
    plt.tight_layout()
    plt.savefig(out_dir / "06_depth_vs_diversity.png", dpi=150, bbox_inches="tight")
    plt.close(); print("[PLOT] 06_depth_vs_diversity.png")


def plot_flagged_heatmap(df_flag, df_meta_summ, out_dir, cov_threshold):
    if df_flag.empty:
        print(f"[INFO] No genera flagged at coverage >= {cov_threshold}"); return

    pivot = df_flag.pivot_table(
        index="sample_id", columns="genus", values="cov", aggfunc="max", fill_value=0
    )
    col_order = (pivot > 0).sum().sort_values(ascending=False).index
    pivot = pivot[col_order]

    # Sort rows by ancestry then max coverage
    if "ancestry" in df_meta_summ.columns:
        anc_map = df_meta_summ.set_index("sample_id")["ancestry"]
        pivot = pivot.copy()
        pivot["_anc"] = pivot.index.map(anc_map).fillna("Unknown")
        pivot = pivot.sort_values("_anc").drop(columns="_anc")

    row_colors = None
    if "ancestry" in df_meta_summ.columns:
        anc_map = df_meta_summ.set_index("sample_id")["ancestry"]
        row_colors = pd.Series(
            [ANCESTRY_COLORS.get(anc_map.get(s, "Unknown"), "#757575") for s in pivot.index],
            index=pivot.index
        )

    fig, ax = plt.subplots(figsize=(max(10, len(pivot.columns) * 0.55), max(6, len(pivot) * 0.18)))
    sns.heatmap(pivot, cmap="YlOrRd", linewidths=0.3, linecolor="#f0f0f0",
                ax=ax, cbar_kws={"label": "Coverage", "shrink": 0.6}, yticklabels=True)
    ax.set_title(f"High-Confidence Genera (coverage ≥ {cov_threshold}) — Coverage per Sample",
                 fontsize=11, fontweight="bold")
    ax.tick_params(axis="y", labelsize=6)
    plt.xticks(rotation=45, ha="right", fontsize=8)

    # Ancestry colour bar on y-axis
    if row_colors is not None:
        for i, sid in enumerate(pivot.index):
            color = row_colors.get(sid, "#757575")
            ax.add_patch(plt.Rectangle((-0.6, i), 0.5, 1, color=color,
                                        transform=ax.get_yaxis_transform(), clip_on=False))
        legend_handles = [mpatches.Patch(color=c, label=a)
                          for a, c in ANCESTRY_COLORS.items() if a != "Unknown"]
        ax.legend(handles=legend_handles, loc="upper right", fontsize=7,
                  title="Ancestry", bbox_to_anchor=(1.18, 1))

    plt.tight_layout()
    plt.savefig(out_dir / "07_flagged_genera_heatmap.png", dpi=150, bbox_inches="tight")
    plt.close(); print("[PLOT] 07_flagged_genera_heatmap.png")


def plot_kmer_matrix(df_kmer, df_meta_summ, out_dir, top_n=40):
    if df_kmer.empty: return

    human_like = [c for c in df_kmer.columns
                  if any(h.lower() in c.lower() for h in ["homo", "human"])]
    df_filt = df_kmer.drop(columns=human_like, errors="ignore")
    top_genera = df_filt.sum(axis=0).nlargest(top_n).index
    df_top = df_filt[top_genera]
    df_log = np.log10(df_top + 1)

    # Cluster rows (only if >2 samples)
    if len(df_log) > 2:
        row_link = linkage(pdist(df_log.fillna(0), metric="braycurtis"), method="average")
        row_order = dendrogram(row_link, no_plot=True)["leaves"]
        df_log = df_log.iloc[row_order]

    fig, ax = plt.subplots(figsize=(max(12, top_n * 0.45), max(8, len(df_log) * 0.14)))
    sns.heatmap(df_log, cmap="Blues", linewidths=0, ax=ax,
                cbar_kws={"label": "log₁₀(k-mers + 1)", "shrink": 0.6}, yticklabels=False)
    ax.set_title(f"Genus K-mer Abundance Matrix — Top {top_n} Genera (log₁₀, Bray-Curtis clustered)",
                 fontsize=11, fontweight="bold")
    ax.set_xlabel("Genus"); ax.set_ylabel(f"Samples (n={len(df_log)})")
    plt.xticks(rotation=45, ha="right", fontsize=8)

    # Ancestry row colour strip
    if "ancestry" in df_meta_summ.columns:
        anc_map = df_meta_summ.set_index("sample_id")["ancestry"]
        for i, sid in enumerate(df_log.index):
            anc = anc_map.get(sid, "Unknown")
            color = ANCESTRY_COLORS.get(anc, "#757575")
            ax.add_patch(plt.Rectangle((-0.8, i), 0.7, 1, color=color,
                                        transform=ax.get_yaxis_transform(), clip_on=False))
        legend_handles = [mpatches.Patch(color=c, label=a)
                          for a, c in ANCESTRY_COLORS.items() if a != "Unknown"]
        ax.legend(handles=legend_handles, loc="upper right", fontsize=7,
                  title="Ancestry", bbox_to_anchor=(1.18, 1))

    plt.tight_layout()
    plt.savefig(out_dir / "08_kmer_matrix_heatmap.png", dpi=150, bbox_inches="tight")
    plt.close(); print("[PLOT] 08_kmer_matrix_heatmap.png")


def plot_flagged_prevalence(df_flag, out_dir, cov_threshold):
    if df_flag.empty: return
    prev = (df_flag.groupby("genus")["sample_id"].nunique()
            .sort_values(ascending=False).head(25).reset_index())
    prev.columns = ["genus", "n_samples"]
    colors = [
        "#FF6F00" if any(c.lower() in g.lower() for c in CONTAMINANTS) else "#0288D1"
        for g in prev["genus"]
    ]
    fig, ax = plt.subplots(figsize=(10, max(5, len(prev) * 0.35)))
    ax.barh(prev["genus"][::-1], prev["n_samples"][::-1], color=colors[::-1])
    ax.set_xlabel("Number of Samples")
    ax.set_title(f"Top 25 Flagged Genera by Prevalence (coverage ≥ {cov_threshold})",
                 fontsize=11, fontweight="bold")
    ax.legend(handles=[mpatches.Patch(facecolor="#0288D1", label="Microbial genus"),
                        mpatches.Patch(facecolor="#FF6F00", label="Known contaminant")],
              loc="lower right", fontsize=9)
    plt.tight_layout()
    plt.savefig(out_dir / "09_flagged_prevalence.png", dpi=150, bbox_inches="tight")
    plt.close(); print("[PLOT] 09_flagged_prevalence.png")


def plot_isup_diversity(df, out_dir):
    """Shannon entropy by ISUP grade — only samples with valid ISUP."""
    sub = df.dropna(subset=["isup"]).copy()
    sub["isup"] = sub["isup"].astype(int)
    isup_vals = sorted(sub["isup"].unique())
    if len(isup_vals) < 2:
        print("[SKIP] Not enough ISUP grades for plot"); return

    fig, ax = plt.subplots(figsize=(9, 5))
    data = [sub.loc[sub["isup"] == g, "shannon"].values for g in isup_vals]
    bp = ax.boxplot(data, patch_artist=True,
                    medianprops=dict(color="black", lw=2))
    for patch in bp["boxes"]:
        patch.set_facecolor("#0288D1"); patch.set_alpha(0.7)
    for i, (g, vals) in enumerate(zip(isup_vals, data)):
        jitter = np.random.default_rng(42).uniform(-0.15, 0.15, len(vals))
        ax.scatter(i + 1 + jitter, vals, color="#01579B", alpha=0.5, s=15, zorder=3)

    valid = [d for d in data if len(d) > 1]
    stat, p = safe_kruskal(*valid)
    if stat is not None:
        ax.set_title(f"Shannon Diversity by ISUP Grade\nKruskal-Wallis p={p:.3f}", fontsize=11)
    else:
        ax.set_title("Shannon Diversity by ISUP Grade")

    ax.set_xticklabels([f"ISUP {g}" for g in isup_vals])
    ax.set_ylabel("Shannon Entropy")
    plt.tight_layout()
    plt.savefig(out_dir / "10_isup_diversity.png", dpi=150, bbox_inches="tight")
    plt.close(); print("[PLOT] 10_isup_diversity.png")


def plot_cohort_summary(df, out_dir):
    """Two-panel: ancestry composition pie + country/region bar."""
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    fig.suptitle("HEROIC1K — Cohort Composition", fontsize=12, fontweight="bold")

    # Pie: ancestry
    anc_counts = df["ancestry"].value_counts()
    colors_pie  = [ANCESTRY_COLORS.get(a, "#757575") for a in anc_counts.index]
    axes[0].pie(anc_counts.values, labels=anc_counts.index, colors=colors_pie,
                autopct="%1.1f%%", startangle=140, pctdistance=0.8)
    axes[0].set_title("Ancestry Distribution")

    # Bar: ethnicity (top 10)
    eth_counts = df["ethnicity"].value_counts().head(10)
    axes[1].barh(eth_counts.index[::-1], eth_counts.values[::-1], color="#5C6BC0")
    axes[1].set_xlabel("Number of Samples")
    axes[1].set_title("Self-Reported Ethnicity (top 10)")

    plt.tight_layout()
    plt.savefig(out_dir / "00_cohort_summary.png", dpi=150, bbox_inches="tight")
    plt.close(); print("[PLOT] 00_cohort_summary.png")


# ── Text summary ───────────────────────────────────────────────────────────────

def write_report(df, df_div, df_flag, df_kmer, cov_threshold, out_dir):
    lines = []
    sep = "=" * 68

    lines += [sep, "HEROIC1K KrakenUniq Evaluation — Summary Report", sep, ""]
    lines.append(f"  Samples analysed:        {len(df)}")
    lines.append(f"  Coverage flag threshold:  >= {cov_threshold}")
    lines.append(f"  Total genera detected:    {df_kmer.shape[1]}")
    lines.append("")

    # Cohort breakdown
    lines.append("── Cohort Composition " + "─" * 39)
    for anc, cnt in df["ancestry"].value_counts().items():
        lines.append(f"  {anc:<25} {cnt:>4} samples")
    lines.append(f"  {'Country — South Africa':<25} {(df['country']=='South Africa').sum():>4} samples")
    lines.append(f"  {'Country — Australia':<25} {(df['country']=='Australia').sum():>4} samples")
    lines.append("")

    lines.append("── Read Depth " + "─" * 47)
    lines.append(f"  Median total reads:      {df['total_reads'].median():>12,.0f}")
    lines.append(f"  Mean   total reads:      {df['total_reads'].mean():>12,.0f}")
    lines.append(f"  Min / Max:               {df['total_reads'].min():>12,.0f} / {df['total_reads'].max():,.0f}")
    lines.append(f"  Samples < 1M reads:      {(df['total_reads'] < 1e6).sum():>4}")
    lines.append("")

    lines.append("── Human Read Fraction " + "─" * 38)
    lines.append(f"  Median human %:          {df['human_pct'].median():>10.2f}%")
    lines.append(f"  Mean   human %:          {df['human_pct'].mean():>10.2f}%")
    lines.append(f"  Std:                     {df['human_pct'].std():>10.2f}%")
    lines.append(f"  Min / Max:               {df['human_pct'].min():>10.2f}% / {df['human_pct'].max():.2f}%")
    lines.append("")

    # Human % by ancestry
    lines.append("  Human % by ancestry:")
    for anc in ["Southern African", "European", "Admixed", "Asian"]:
        sub = df[df["ancestry"] == anc]["human_pct"]
        if sub.empty: continue
        lines.append(f"    {anc:<20}  median={sub.median():.2f}%  mean={sub.mean():.2f}%  n={len(sub)}")
    lines.append("")

    lines.append("── Unclassified Fraction " + "─" * 36)
    lines.append(f"  Median unclassified %:   {df['unclassified_pct'].median():>10.2f}%")
    lines.append(f"  Mean   unclassified %:   {df['unclassified_pct'].mean():>10.2f}%")
    lines.append("")

    lines.append("── Alpha Diversity (Genus Level) " + "─" * 28)
    lines.append(f"  Median observed genera:  {df_div['observed_genera'].median():>8.0f}")
    lines.append(f"  Median Shannon entropy:  {df_div['shannon'].median():>10.3f}")
    lines.append(f"  Median Simpson index:    {df_div['simpson'].median():>10.3f}")
    r, p = spearmanr(df["total_reads"], df_div["shannon"])
    lines.append(f"  Spearman r (depth~shannon): {r:.3f}  p={p:.3e}")
    lines.append("")

    lines.append(f"── Flagged Genera (coverage >= {cov_threshold}) " + "─" * 24)
    if df_flag.empty:
        lines.append("  None detected.")
    else:
        lines.append(f"  Samples with >= 1 flagged genus:  {df_flag['sample_id'].nunique()}")
        lines.append(f"  Unique flagged genera:            {df_flag['genus'].nunique()}")
        lines.append("")
        lines.append("  Top flagged genera by prevalence:")
        prev = df_flag.groupby("genus").agg(
            n_samples=("sample_id","nunique"),
            median_cov=("cov","median"),
            median_kmers=("kmers","median"),
        ).sort_values("n_samples", ascending=False).head(20)
        lines.append(f"  {'Genus':<42} {'Samples':>8} {'Med.Cov':>10} {'Med.Kmers':>12}")
        lines.append("  " + "-"*74)
        for genus, row in prev.iterrows():
            tag = " [CONTAM]" if any(c.lower() in genus.lower() for c in CONTAMINANTS) else ""
            lines.append(
                f"  {(genus+tag):<42} {int(row['n_samples']):>8} "
                f"{row['median_cov']:>10.5f} {int(row['median_kmers']):>12,}"
            )
    lines.append("")

    lines.append("── Outlier / QC Notes " + "─" * 39)
    hi = df[df["human_pct"] > 95]
    lo = df[df["total_reads"] < 1_000_000]
    lines.append(f"  Samples >95% human:   {len(hi)}")
    lines.append(f"  Low-depth (<1M):      {len(lo)}")
    if not lo.empty:
        for _, r in lo.iterrows():
            anc = r.get("ancestry", "?")
            lines.append(f"    {r['sample_id']:<20} {r['total_reads']:>10,}  [{anc}]")

    lines += ["", sep]
    text = "\n".join(lines)
    print(text)
    with open(out_dir / "summary_report.txt", "w") as fh:
        fh.write(text + "\n")
    print("\n[OUT] summary_report.txt")


# ── Main ───────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--report_dir",    required=True)
    ap.add_argument("--metadata",      default=None,
                    help="Path to GEOETHNIC_MASTER_BATCH1.csv")
    ap.add_argument("--output_dir",    default="./heroic1k_results")
    ap.add_argument("--pattern",       default="*.krakenuniq.report")
    ap.add_argument("--cov_threshold", type=float, default=0.005)
    ap.add_argument("--top_genera",    type=int,   default=40)
    args = ap.parse_args()

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    print(f"[INFO] Output:   {out_dir.resolve()}")
    print(f"[INFO] Cov flag: >= {args.cov_threshold}")

    # Load reports
    df_summ, df_div, df_flag, df_kmer = load_all(
        args.report_dir, args.pattern, args.cov_threshold)

    # Load + merge metadata
    if args.metadata:
        meta = load_metadata(args.metadata)
        df_merged = merge_metadata(df_summ, meta)
        df_div_merged = df_div.set_index("sample_id").join(
            df_merged.set_index("sample_id")[["ancestry","country","isup","age","psa","total_reads"]],
            how="left"
        ).reset_index()
    else:
        df_merged      = df_summ.copy()
        df_merged["ancestry"] = "Unknown"
        df_merged["country"]  = "Unknown"
        df_div_merged = df_div.merge(df_summ[["sample_id","total_reads"]], on="sample_id", how="left")
        df_div_merged["ancestry"] = "Unknown"
        print("[WARN] No metadata provided — ancestry plots will be uninformative")

    # Save tables
    df_merged.to_csv(out_dir / "per_sample_summary.csv", index=False)
    df_div_merged.to_csv(out_dir / "alpha_diversity.csv", index=False)
    df_kmer.to_csv(out_dir / "genus_kmer_matrix.csv")
    if not df_flag.empty:
        df_flag.to_csv(out_dir / "flagged_genera.csv", index=False)
        df_flag.pivot_table(index="sample_id", columns="genus",
                            values="cov", aggfunc="max", fill_value=0
                            ).to_csv(out_dir / "flagged_cov_pivot.csv")
        df_flag.pivot_table(index="sample_id", columns="genus",
                            values="kmers", aggfunc="max", fill_value=0
                            ).to_csv(out_dir / "flagged_kmer_pivot.csv")
    print("[OUT] CSVs written")

    # Plots
    plot_cohort_summary(df_merged, out_dir)
    plot_read_overview(df_merged, out_dir)
    plot_composition_bar(df_merged, out_dir)
    plot_ancestry_read_stats(df_merged, out_dir)
    plot_alpha_diversity(df_div_merged, out_dir)
    plot_ancestry_diversity(df_div_merged, out_dir)
    plot_depth_diversity_scatter(df_div_merged, out_dir)
    plot_flagged_heatmap(df_flag, df_merged, out_dir, args.cov_threshold)
    plot_kmer_matrix(df_kmer, df_merged, out_dir, top_n=args.top_genera)
    plot_flagged_prevalence(df_flag, out_dir, args.cov_threshold)
    plot_isup_diversity(df_div_merged, out_dir)

    write_report(df_merged, df_div_merged, df_flag, df_kmer, args.cov_threshold, out_dir)

    print(f"\n[DONE] Results in: {out_dir.resolve()}")


if __name__ == "__main__":
    main()
