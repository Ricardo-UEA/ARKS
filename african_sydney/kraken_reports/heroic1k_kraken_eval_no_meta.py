#!/usr/bin/env python3
"""
HEROIC1K KrakenUniq Report Evaluation Pipeline
================================================
Parses KrakenUniq reports from the HEROIC1K African prostate cancer cohort.
Focuses on genus-level analysis with:
  - Per-sample summary (read depth, human fraction, unclassified rate)
  - Coverage-flagged taxa: genus-level hits with coverage >= threshold (default 0.005)
  - Genus k-mer abundance matrix (all samples x all genera)
  - Alpha diversity metrics
  - Publication-quality plots

KrakenUniq report columns (tab-separated):
    %  |  reads  |  taxReads  |  kmers  |  dup  |  cov  |  taxID  |  rank  |  taxName

Usage:
    python heroic1k_kraken_eval.py \
        --report_dir /gpfs/.../kraken_reports/no_depletion_reports \
        --output_dir ./heroic1k_results \
        --cov_threshold 0.005
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
import seaborn as sns
from matplotlib.patches import Patch

# ── Constants ─────────────────────────────────────────────────────────────────

HUMAN_NAMES  = {"Homo sapiens", "Homo"}
HUMAN_TAXIDS = {"9606", "9605"}

CONTAMINANTS = {"Buchnera", "Lymphocryptovirus", "Epstein-Barr virus"}

PALETTE = {
    "human":     "#E53935",
    "nonhuman":  "#43A047",
    "unclass":   "#9E9E9E",
    "flagged":   "#FF6F00",
    "diversity": ["#7B1FA2", "#0288D1", "#F4511E"],
}

# ── Parser ─────────────────────────────────────────────────────────────────────

def parse_report(filepath):
    """
    Parse a single KrakenUniq report.
    Returns:
        genera        : list of genus-level row dicts
        total_reads   : int (classified + unclassified)
        human_reads   : int
        unclass_reads : int
        classified    : int
    """
    genera = []
    unclass_reads = 0
    root_reads    = 0
    human_reads   = 0

    with open(filepath) as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            if parts[0].strip() == "%":
                continue  # header row

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
                    "reads_clade":  reads_clade,
                    "reads_direct": reads_direct,
                    "kmers":        kmers,
                    "dup":          dup,
                    "cov":          cov,
                    "taxid":        taxid,
                    "name":         name,
                })

    total_reads = unclass_reads + root_reads
    return genera, total_reads, human_reads, unclass_reads, root_reads


# ── Per-sample summary ─────────────────────────────────────────────────────────

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


# ── Diversity helpers ──────────────────────────────────────────────────────────

def shannon(counts):
    a = np.array(counts, dtype=float)
    a = a[a > 0]
    if a.sum() == 0:
        return 0.0
    p = a / a.sum()
    return float(-np.sum(p * np.log(p)))

def simpson(counts):
    a = np.array(counts, dtype=float)
    n = a.sum()
    if n < 2:
        return 0.0
    return float(1 - np.sum(a * (a - 1)) / (n * (n - 1)))


# ── Load all reports ───────────────────────────────────────────────────────────

def load_all(report_dir, pattern, cov_threshold):
    report_dir = Path(report_dir)
    files = sorted(report_dir.glob(pattern))
    if not files:
        files = sorted(report_dir.glob("*.report"))
    if not files:
        sys.exit(f"[ERROR] No report files found in {report_dir}")

    print(f"[INFO] Found {len(files)} report files")

    summaries   = []
    diversity   = []
    flagged     = []
    kmer_matrix = defaultdict(dict)   # genus_name -> {sample_id: kmers}

    for f in files:
        sid = f.name.split(".krakenuniq")[0].split(".report")[0]
        genera, total, human, unclass, classified = parse_report(f)

        summaries.append(summarise_sample(sid, genera, total, human, unclass, classified))

        kmer_counts = [g["kmers"] for g in genera if g["kmers"] > 0]
        diversity.append({
            "sample_id":       sid,
            "observed_genera": len(genera),
            "shannon":         shannon(kmer_counts),
            "simpson":         simpson(kmer_counts),
        })

        for g in genera:
            name = g["name"]
            kmer_matrix[name][sid] = g["kmers"]

            if g["cov"] >= cov_threshold:
                flagged.append({
                    "sample_id":      sid,
                    "genus":          name,
                    "reads_clade":    g["reads_clade"],
                    "kmers":          g["kmers"],
                    "dup":            g["dup"],
                    "cov":            g["cov"],
                    "is_contaminant": any(c.lower() in name.lower() for c in CONTAMINANTS),
                })

    df_summ = pd.DataFrame(summaries)
    df_div  = pd.DataFrame(diversity)
    df_flag = pd.DataFrame(flagged) if flagged else pd.DataFrame(
        columns=["sample_id","genus","reads_clade","kmers","dup","cov","is_contaminant"])

    all_samples = [s["sample_id"] for s in summaries]
    df_kmer = pd.DataFrame(kmer_matrix, index=all_samples).fillna(0).astype(int)

    return df_summ, df_div, df_flag, df_kmer


# ── Plots ──────────────────────────────────────────────────────────────────────

def plot_read_overview(df, out_dir):
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    fig.suptitle("HEROIC1K — Read-Level Overview (No-Depletion Baseline)",
                 fontsize=13, fontweight="bold")

    configs = [
        ("human_pct",        "Human Reads (%)",       PALETTE["human"]),
        ("unclassified_pct", "Unclassified Reads (%)", PALETTE["unclass"]),
        ("total_reads",      "Total Reads (M)",        "#0288D1"),
    ]
    for ax, (col, label, color) in zip(axes, configs):
        vals = df[col] / (1e6 if col == "total_reads" else 1)
        ax.hist(vals, bins=30, color=color, edgecolor="white", linewidth=0.4)
        med = vals.median()
        ax.axvline(med, color="black", ls="--", lw=1.2, label=f"Median: {med:.2f}")
        ax.set_xlabel(label)
        ax.set_ylabel("Samples")
        ax.set_title(label)
        ax.legend(fontsize=9)

    plt.tight_layout()
    plt.savefig(out_dir / "01_read_overview.png", dpi=150, bbox_inches="tight")
    plt.close()
    print("[PLOT] 01_read_overview.png")


def plot_composition_bar(df, out_dir):
    df_s  = df.sort_values("human_pct", ascending=False).reset_index(drop=True)
    total = df_s["total_reads"].replace(0, 1)
    h  = df_s["human_reads"] / total * 100
    nc = df_s["nonhuman_classified"] / total * 100
    u  = df_s["unclassified_reads"] / total * 100

    fig, ax = plt.subplots(figsize=(max(14, len(df_s) * 0.1), 5))
    x = np.arange(len(df_s))
    ax.bar(x, h,            label="Human",                   color=PALETTE["human"],    width=1.0)
    ax.bar(x, nc, bottom=h, label="Non-human classified",    color=PALETTE["nonhuman"], width=1.0)
    ax.bar(x, u,  bottom=h+nc, label="Unclassified",         color=PALETTE["unclass"],  width=1.0)
    ax.set_xlim(-0.5, len(df_s) - 0.5)
    ax.set_ylim(0, 100)
    ax.set_xlabel(f"Samples (n={len(df_s)}, sorted by human %)")
    ax.set_ylabel("Read Composition (%)")
    ax.set_title("HEROIC1K — Read Composition per Sample")
    ax.set_xticks([])
    ax.legend(loc="upper right")
    plt.tight_layout()
    plt.savefig(out_dir / "02_read_composition.png", dpi=150, bbox_inches="tight")
    plt.close()
    print("[PLOT] 02_read_composition.png")


def plot_alpha_diversity(df_div, out_dir):
    fig, axes = plt.subplots(1, 3, figsize=(15, 4))
    fig.suptitle("HEROIC1K — Genus-Level Alpha Diversity", fontsize=12, fontweight="bold")

    configs = [
        ("observed_genera", "Observed Genera"),
        ("shannon",         "Shannon Entropy (k-mers)"),
        ("simpson",         "Simpson Index (k-mers)"),
    ]
    for ax, (col, label), color in zip(axes, configs, PALETTE["diversity"]):
        ax.hist(df_div[col], bins=25, color=color, edgecolor="white", linewidth=0.4)
        med = df_div[col].median()
        ax.axvline(med, color="black", ls="--", lw=1.2, label=f"Median: {med:.2f}")
        ax.set_xlabel(label)
        ax.set_ylabel("Samples")
        ax.set_title(label)
        ax.legend(fontsize=9)

    plt.tight_layout()
    plt.savefig(out_dir / "03_alpha_diversity.png", dpi=150, bbox_inches="tight")
    plt.close()
    print("[PLOT] 03_alpha_diversity.png")


def plot_flagged_heatmap(df_flag, out_dir, cov_threshold):
    if df_flag.empty:
        print(f"[INFO] No genera flagged at coverage >= {cov_threshold} — skipping heatmap")
        return

    pivot = df_flag.pivot_table(
        index="sample_id", columns="genus", values="cov",
        aggfunc="max", fill_value=0
    )
    col_order = (pivot > 0).sum().sort_values(ascending=False).index
    pivot = pivot[col_order]

    fig, ax = plt.subplots(figsize=(max(10, len(pivot.columns) * 0.55), max(6, len(pivot) * 0.18)))
    sns.heatmap(
        pivot, cmap="YlOrRd", linewidths=0.3, linecolor="#f0f0f0",
        ax=ax, cbar_kws={"label": "Coverage", "shrink": 0.6}, yticklabels=True,
    )
    ax.set_title(
        f"High-Confidence Genera (coverage ≥ {cov_threshold}) — Coverage per Sample",
        fontsize=11, fontweight="bold"
    )
    ax.set_xlabel("Genus")
    ax.set_ylabel("Sample")
    ax.tick_params(axis="y", labelsize=6)
    plt.xticks(rotation=45, ha="right", fontsize=8)
    plt.tight_layout()
    plt.savefig(out_dir / "04_flagged_genera_heatmap.png", dpi=150, bbox_inches="tight")
    plt.close()
    print("[PLOT] 04_flagged_genera_heatmap.png")


def plot_kmer_matrix(df_kmer, out_dir, top_n=40):
    if df_kmer.empty:
        print("[WARN] K-mer matrix empty — skipping")
        return

    human_like = [c for c in df_kmer.columns
                  if any(h.lower() in c.lower() for h in ["homo", "human"])]
    df_filt = df_kmer.drop(columns=human_like, errors="ignore")

    top_genera = df_filt.sum(axis=0).nlargest(top_n).index
    df_top = df_filt[top_genera]
    df_log = np.log10(df_top + 1)

    fig, ax = plt.subplots(figsize=(max(12, top_n * 0.45), max(8, len(df_log) * 0.16)))
    sns.heatmap(
        df_log, cmap="Blues", linewidths=0, ax=ax,
        cbar_kws={"label": "log₁₀(k-mers + 1)", "shrink": 0.6},
        yticklabels=False,
    )
    ax.set_title(
        f"Genus K-mer Abundance Matrix — Top {top_n} Genera (log₁₀ scale)",
        fontsize=11, fontweight="bold"
    )
    ax.set_xlabel("Genus")
    ax.set_ylabel(f"Samples (n={len(df_log)})")
    plt.xticks(rotation=45, ha="right", fontsize=8)
    plt.tight_layout()
    plt.savefig(out_dir / "05_kmer_matrix_heatmap.png", dpi=150, bbox_inches="tight")
    plt.close()
    print("[PLOT] 05_kmer_matrix_heatmap.png")


def plot_flagged_prevalence(df_flag, out_dir, cov_threshold):
    if df_flag.empty:
        return

    prevalence = (
        df_flag.groupby("genus")["sample_id"].nunique()
        .sort_values(ascending=False).head(25).reset_index()
    )
    prevalence.columns = ["genus", "n_samples"]

    colors = [
        PALETTE["flagged"] if any(c.lower() in g.lower() for c in CONTAMINANTS) else "#0288D1"
        for g in prevalence["genus"]
    ]

    fig, ax = plt.subplots(figsize=(10, max(5, len(prevalence) * 0.35)))
    ax.barh(prevalence["genus"][::-1], prevalence["n_samples"][::-1], color=colors[::-1])
    ax.set_xlabel("Number of Samples")
    ax.set_title(
        f"Top 25 Flagged Genera by Prevalence (coverage ≥ {cov_threshold})",
        fontsize=11, fontweight="bold"
    )
    legend_elements = [
        Patch(facecolor="#0288D1",        label="Microbial genus"),
        Patch(facecolor=PALETTE["flagged"], label="Known contaminant"),
    ]
    ax.legend(handles=legend_elements, loc="lower right", fontsize=9)
    plt.tight_layout()
    plt.savefig(out_dir / "06_flagged_prevalence.png", dpi=150, bbox_inches="tight")
    plt.close()
    print("[PLOT] 06_flagged_prevalence.png")


# ── Text report ────────────────────────────────────────────────────────────────

def write_report(df_summ, df_div, df_flag, df_kmer, cov_threshold, out_dir):
    lines = []
    sep = "=" * 68

    lines += [sep, "HEROIC1K KrakenUniq Evaluation — Summary Report", sep, ""]
    lines.append(f"  Samples analysed:        {len(df_summ)}")
    lines.append(f"  Coverage flag threshold:  >= {cov_threshold}")
    lines.append(f"  Total genera detected:    {df_kmer.shape[1]}")
    lines.append("")

    lines.append("── Read Depth " + "─" * 47)
    lines.append(f"  Median total reads:      {df_summ['total_reads'].median():>12,.0f}")
    lines.append(f"  Mean   total reads:      {df_summ['total_reads'].mean():>12,.0f}")
    lines.append(f"  Min / Max:               {df_summ['total_reads'].min():>12,.0f} / {df_summ['total_reads'].max():,.0f}")
    lines.append("")

    lines.append("── Human Read Fraction " + "─" * 38)
    lines.append(f"  Median human %:          {df_summ['human_pct'].median():>10.2f}%")
    lines.append(f"  Mean   human %:          {df_summ['human_pct'].mean():>10.2f}%")
    lines.append(f"  Std:                     {df_summ['human_pct'].std():>10.2f}%")
    lines.append(f"  Min / Max:               {df_summ['human_pct'].min():>10.2f}% / {df_summ['human_pct'].max():.2f}%")
    lines.append("")

    lines.append("── Unclassified Fraction " + "─" * 36)
    lines.append(f"  Median unclassified %:   {df_summ['unclassified_pct'].median():>10.2f}%")
    lines.append(f"  Mean   unclassified %:   {df_summ['unclassified_pct'].mean():>10.2f}%")
    lines.append("")

    lines.append("── Alpha Diversity (Genus Level) " + "─" * 28)
    lines.append(f"  Median observed genera:  {df_div['observed_genera'].median():>8.0f}")
    lines.append(f"  Median Shannon entropy:  {df_div['shannon'].median():>10.3f}")
    lines.append(f"  Median Simpson index:    {df_div['simpson'].median():>10.3f}")
    lines.append("")

    lines.append(f"── Flagged Genera (coverage >= {cov_threshold}) " + "─" * 26)
    if df_flag.empty:
        lines.append("  None detected.")
    else:
        n_samp  = df_flag["sample_id"].nunique()
        n_gen   = df_flag["genus"].nunique()
        lines.append(f"  Samples with >= 1 flagged genus:  {n_samp}")
        lines.append(f"  Unique flagged genera:            {n_gen}")
        lines.append("")
        lines.append("  Top flagged genera by prevalence:")
        prev = df_flag.groupby("genus").agg(
            n_samples    =("sample_id", "nunique"),
            median_cov   =("cov",       "median"),
            median_kmers =("kmers",     "median"),
        ).sort_values("n_samples", ascending=False).head(20)
        lines.append(f"  {'Genus':<42} {'Samples':>8} {'Med.Cov':>10} {'Med.Kmers':>12}")
        lines.append("  " + "-" * 74)
        for genus, row in prev.iterrows():
            tag = " [CONTAM]" if any(c.lower() in genus.lower() for c in CONTAMINANTS) else ""
            lines.append(
                f"  {(genus+tag):<42} {int(row['n_samples']):>8} "
                f"{row['median_cov']:>10.5f} {int(row['median_kmers']):>12,}"
            )
    lines.append("")

    lines.append("── Per-Sample Flagged Notes " + "─" * 33)
    if df_flag.empty:
        lines.append("  None.")
    else:
        per_sample = df_flag.groupby("sample_id").apply(
            lambda x: pd.Series({
                "n_flagged": x["genus"].nunique(),
                "max_cov":   x["cov"].max(),
                "top_genus": x.loc[x["cov"].idxmax(), "genus"],
            })
        ).sort_values("max_cov", ascending=False)
        lines.append(f"  {'Sample':<20} {'#Flagged':>10} {'Max Cov':>10}  Top Genus")
        lines.append("  " + "-" * 68)
        for sid, row in per_sample.iterrows():
            lines.append(
                f"  {sid:<20} {int(row['n_flagged']):>10} "
                f"{row['max_cov']:>10.5f}  {row['top_genus']}"
            )
    lines.append("")

    lines.append("── Outlier Samples " + "─" * 41)
    hi = df_summ[df_summ["human_pct"] > 95]
    lo = df_summ[df_summ["total_reads"] < 1_000_000]
    lines.append("  Samples >95% human reads:" if not hi.empty else "  No samples >95% human reads.")
    for _, r in hi.iterrows():
        lines.append(f"    {r['sample_id']}: {r['human_pct']:.1f}%")
    lines.append("  Low-depth samples (<1M reads):" if not lo.empty else "  No low-depth samples (<1M).")
    for _, r in lo.iterrows():
        lines.append(f"    {r['sample_id']}: {r['total_reads']:,}")

    lines += ["", sep]
    text = "\n".join(lines)
    print(text)
    with open(out_dir / "summary_report.txt", "w") as fh:
        fh.write(text + "\n")
    print("\n[OUT] summary_report.txt")


# ── Main ───────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(description="HEROIC1K KrakenUniq genus-level evaluation")
    ap.add_argument("--report_dir",    required=True)
    ap.add_argument("--output_dir",    default="./heroic1k_results")
    ap.add_argument("--pattern",       default="*.krakenuniq.report",
                    help="Glob pattern for report files")
    ap.add_argument("--cov_threshold", type=float, default=0.005,
                    help="Coverage threshold to flag high-confidence genera (default: 0.005)")
    ap.add_argument("--top_genera",    type=int, default=40,
                    help="Top N genera for k-mer matrix heatmap (default: 40)")
    args = ap.parse_args()

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    print(f"[INFO] Output:    {out_dir.resolve()}")
    print(f"[INFO] Cov flag:  >= {args.cov_threshold}")

    df_summ, df_div, df_flag, df_kmer = load_all(
        args.report_dir, args.pattern, args.cov_threshold
    )

    # Save tables
    df_summ.to_csv(out_dir / "per_sample_summary.csv", index=False)
    df_div.to_csv(out_dir  / "alpha_diversity.csv", index=False)
    df_kmer.to_csv(out_dir / "genus_kmer_matrix.csv")
    if not df_flag.empty:
        df_flag.to_csv(out_dir / "flagged_genera.csv", index=False)
        df_flag.pivot_table(
            index="sample_id", columns="genus", values="cov", aggfunc="max", fill_value=0
        ).to_csv(out_dir / "flagged_cov_pivot.csv")
        df_flag.pivot_table(
            index="sample_id", columns="genus", values="kmers", aggfunc="max", fill_value=0
        ).to_csv(out_dir / "flagged_kmer_pivot.csv")
    print("[OUT] CSVs written")

    # Plots
    plot_read_overview(df_summ, out_dir)
    plot_composition_bar(df_summ, out_dir)
    plot_alpha_diversity(df_div, out_dir)
    plot_flagged_heatmap(df_flag, out_dir, args.cov_threshold)
    plot_kmer_matrix(df_kmer, out_dir, top_n=args.top_genera)
    plot_flagged_prevalence(df_flag, out_dir, args.cov_threshold)

    write_report(df_summ, df_div, df_flag, df_kmer, args.cov_threshold, out_dir)

    print(f"\n[DONE] Results in: {out_dir.resolve()}")


if __name__ == "__main__":
    main()
