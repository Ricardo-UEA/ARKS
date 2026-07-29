#!/usr/bin/env python3
"""
HEROIC1K — No-Depletion vs Post-Depletion Comparison
======================================================
Reads pre-computed CSVs from both result directories and generates
direct before/after comparison figures.

Usage:
    python heroic1k_compare_t2t.py \
        --no_dep_dir  /gpfs/.../heroic1k_results_no_depletion/meta_data_analysis_included \
        --dep_dir     /gpfs/.../heroic1k_results_depletion \
        --output_dir  /gpfs/home/hce24xau/scratch/gen_kmers/data/african_sydney/kraken_reports/heroic1k_results/heroic1k_comparison
"""

import argparse
from pathlib import Path

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
import seaborn as sns
from scipy.stats import mannwhitneyu, wilcoxon
from statsmodels.stats.multitest import multipletests

# ── Palette ───────────────────────────────────────────────────────────────────

ANCESTRY_COLORS = {
    "Southern African": "#E53935",
    "European":         "#1565C0",
    "Admixed":          "#6A1B9A",
    "Asian":            "#2E7D32",
    "Unknown":          "#757575",
}
ANC_ORDER = ["Southern African", "European", "Admixed", "Asian", "Unknown"]

C_PRE  = "#607D8B"   # blue-grey  — no depletion
C_POST = "#00897B"   # teal       — post depletion

# ── Loaders ───────────────────────────────────────────────────────────────────

def load_dir(result_dir):
    d = Path(result_dir)
    summ = pd.read_csv(d / "per_sample_summary.csv")
    div  = pd.read_csv(d / "alpha_diversity.csv")
    # Merge diversity into summary for convenience
    key_cols = ["sample_id", "observed_genera", "shannon", "simpson"]
    merged = summ.merge(div[key_cols], on="sample_id", how="left")
    return merged


def align(pre, post):
    """Keep only samples present in both, aligned by sample_id."""
    common = set(pre["sample_id"]) & set(post["sample_id"])
    print(f"[INFO] Samples in both conditions: {len(common)}")
    pre  = pre[pre["sample_id"].isin(common)].set_index("sample_id").sort_index()
    post = post[post["sample_id"].isin(common)].set_index("sample_id").sort_index()
    return pre, post


# ── Plot helpers ──────────────────────────────────────────────────────────────

def ancestry_legend():
    return [mpatches.Patch(color=c, label=a)
            for a, c in ANCESTRY_COLORS.items() if a != "Unknown"]


def mw_pval(a, b):
    """Two-sided Mann-Whitney U, returns p-value string."""
    try:
        _, p = mannwhitneyu(a, b, alternative="two-sided")
        return p
    except Exception:
        return float("nan")


def pairwise_delta_posthoc(df, metric, anc_order):
    """
    All pairwise Mann-Whitney U tests on delta values between ancestry groups.
    Returns a DataFrame with raw and BH-corrected p-values.
    """
    from itertools import combinations
    anc_present = [a for a in anc_order if a in df["ancestry"].values]
    pairs, stats, pvals = [], [], []
    for a1, a2 in combinations(anc_present, 2):
        v1 = df.loc[df["ancestry"] == a1, metric].values
        v2 = df.loc[df["ancestry"] == a2, metric].values
        if len(v1) < 3 or len(v2) < 3:
            continue
        stat, p = mannwhitneyu(v1, v2, alternative="two-sided")
        pairs.append((a1, a2))
        stats.append(stat)
        pvals.append(p)

    if not pvals:
        return pd.DataFrame()

    reject, p_adj, _, _ = multipletests(pvals, method="fdr_bh")
    rows = []
    for (a1, a2), stat, p_raw, p_cor, rej in zip(pairs, stats, pvals, p_adj, reject):
        v1 = df.loc[df["ancestry"] == a1, metric].values
        v2 = df.loc[df["ancestry"] == a2, metric].values
        rows.append({
            "group1":      a1,
            "group2":      a2,
            "n1":          len(v1),
            "n2":          len(v2),
            "median1":     np.median(v1),
            "median2":     np.median(v2),
            "delta_medians": np.median(v2) - np.median(v1),
            "U_stat":      stat,
            "p_raw":       p_raw,
            "p_adj_BH":    p_cor,
            "significant": rej,
        })
    return pd.DataFrame(rows)


def wilcoxon_pval(a, b):
    """Paired Wilcoxon signed-rank, returns p-value."""
    diff = np.array(a) - np.array(b)
    diff = diff[diff != 0]
    if len(diff) < 10:
        return float("nan")
    _, p = wilcoxon(diff)
    return p


# ── Plot 1: Human read % waterfall before/after ───────────────────────────────

def plot_human_waterfall(pre, post, out_dir, label="PanHuman"):
    df = pd.DataFrame({
        "sample_id": pre.index,
        "pre":       pre["human_pct"].values,
        "post":      post["human_pct"].values,
        "ancestry":  pre["ancestry"].fillna("Unknown").values,
    }).sort_values("pre", ascending=False).reset_index(drop=True)

    fig, axes = plt.subplots(2, 1, figsize=(16, 8), sharex=True)
    fig.suptitle(f"Human Read Fraction — Before vs After {label} Depletion",
                 fontsize=13, fontweight="bold")

    x = np.arange(len(df))
    colors = [ANCESTRY_COLORS.get(a, "#757575") for a in df["ancestry"]]

    axes[0].bar(x, df["pre"],  color=colors, width=1.0)
    axes[0].set_ylabel("Human Reads (%) — Pre-depletion")
    axes[0].set_ylim(0, df["pre"].max() * 1.1)

    axes[1].bar(x, df["post"], color=colors, width=1.0)
    axes[1].set_ylabel("Human Reads (%) — Post-depletion")
    axes[1].set_ylim(0, max(df["post"].max() * 1.1, 0.01))
    axes[1].set_xlabel(f"Samples (n={len(df)}, sorted by pre-depletion human %)")
    axes[1].set_xticks([])

    # Shared ancestry legend
    axes[0].legend(handles=ancestry_legend(), loc="upper right", fontsize=8,
                   title="Ancestry")

    # Depletion efficiency annotation
    eff = (df["pre"] - df["post"]) / df["pre"].replace(0, np.nan) * 100
    axes[0].set_title(f"Pre-depletion  |  Median: {df['pre'].median():.2f}%", fontsize=10)
    axes[1].set_title(f"Post-depletion |  Median: {df['post'].median():.4f}%  "
                      f"|  Median depletion efficiency: {eff.median():.1f}%", fontsize=10)

    plt.tight_layout()
    plt.savefig(out_dir / "C01_human_waterfall.png", dpi=150, bbox_inches="tight")
    plt.close(); print("[PLOT] C01_human_waterfall.png")


# ── Plot 2: Shannon before/after — paired dots by ancestry ───────────────────

def plot_shannon_paired(pre, post, out_dir, label="PanHuman"):
    anc_present = [a for a in ANC_ORDER if a in pre["ancestry"].values]

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    fig.suptitle(f"Alpha Diversity — Before vs After {label} Depletion",
                 fontsize=13, fontweight="bold")

    for ax, metric, label in zip(
        axes,
        ["shannon", "observed_genera"],
        ["Shannon Entropy", "Observed Genera"],
    ):
        for anc in anc_present:
            mask = pre["ancestry"] == anc
            pre_vals  = pre.loc[mask, metric].values
            post_vals = post.loc[mask, metric].values
            color = ANCESTRY_COLORS.get(anc, "#757575")

            # Draw connecting lines
            for pv, dv in zip(pre_vals, post_vals):
                ax.plot([0, 1], [pv, dv], color=color, alpha=0.15, lw=0.8)

            # Draw group medians as thick lines
            ax.plot([0, 1], [np.median(pre_vals), np.median(post_vals)],
                    color=color, lw=3, zorder=5,
                    label=f"{anc} (n={mask.sum()})")

            ax.scatter([0]*len(pre_vals),  pre_vals,  color=color, alpha=0.4, s=18, zorder=4)
            ax.scatter([1]*len(post_vals), post_vals, color=color, alpha=0.4, s=18, zorder=4)

        # Wilcoxon p across all samples
        p = wilcoxon_pval(pre[metric].values, post[metric].values)
        ax.set_title(f"{label}\nPaired Wilcoxon p={p:.3e}", fontsize=10)
        ax.set_xticks([0, 1])
        ax.set_xticklabels(["Pre-depletion", "Post-depletion"], fontsize=11)
        ax.set_ylabel(label)
        ax.set_xlim(-0.3, 1.3)
        ax.legend(fontsize=7, loc="upper left")

    plt.tight_layout()
    plt.savefig(out_dir / "C02_diversity_paired.png", dpi=150, bbox_inches="tight")
    plt.close(); print("[PLOT] C02_diversity_paired.png")


# ── Plot 3: Side-by-side boxplots per ancestry ────────────────────────────────

def plot_ancestry_sidebyside(pre, post, out_dir, label="PanHuman"):
    anc_present = [a for a in ANC_ORDER if a in pre["ancestry"].values]
    metrics = [
        ("human_pct",        "Human Reads (%)"),
        ("unclassified_pct", "Unclassified Reads (%)"),
        ("shannon",          "Shannon Entropy"),
        ("observed_genera",  "Observed Genera"),
    ]

    fig, axes = plt.subplots(2, 2, figsize=(16, 11))
    fig.suptitle("Pre vs Post-Depletion Metrics by Ancestry", fontsize=13, fontweight="bold")

    for ax, (col, label) in zip(axes.flat, metrics):
        x_positions = []
        x_labels    = []
        tick_pos    = []

        for i, anc in enumerate(anc_present):
            base = i * 3
            pre_vals  = pre.loc[pre["ancestry"] == anc, col].dropna().values
            post_vals = post.loc[post["ancestry"] == anc, col].dropna().values
            color = ANCESTRY_COLORS.get(anc, "#757575")

            bp1 = ax.boxplot(pre_vals,  positions=[base],   widths=0.7,
                             patch_artist=True, notch=False,
                             medianprops=dict(color="black", lw=2))
            bp2 = ax.boxplot(post_vals, positions=[base+1], widths=0.7,
                             patch_artist=True, notch=False,
                             medianprops=dict(color="white", lw=2))

            bp1["boxes"][0].set_facecolor(color); bp1["boxes"][0].set_alpha(0.85)
            bp2["boxes"][0].set_facecolor(color); bp2["boxes"][0].set_alpha(0.4)

            tick_pos.append(base + 0.5)
            x_labels.append(anc.replace(" ", "\n"))

            # Mann-Whitney p between pre and post within ancestry
            if len(pre_vals) > 1 and len(post_vals) > 1:
                p = mw_pval(pre_vals, post_vals)
                ymax = max(np.percentile(pre_vals, 95), np.percentile(post_vals, 95))
                sig = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"
                ax.text(base + 0.5, ymax * 1.05, sig, ha="center", fontsize=9, color=color)

        ax.set_xticks(tick_pos)
        ax.set_xticklabels(x_labels, fontsize=8)
        ax.set_title(label, fontsize=10)
        ax.set_ylabel(label)

    # Shared legend
    pre_patch  = mpatches.Patch(facecolor="#607D8B", alpha=0.85, label="Pre-depletion (solid)")
    post_patch = mpatches.Patch(facecolor="#607D8B", alpha=0.4,  label="Post-depletion (faded)")
    fig.legend(handles=[pre_patch, post_patch], loc="lower center",
               ncol=2, fontsize=9, bbox_to_anchor=(0.5, -0.01))

    plt.tight_layout()
    plt.savefig(out_dir / "C03_ancestry_sidebyside.png", dpi=150, bbox_inches="tight")
    plt.close(); print("[PLOT] C03_ancestry_sidebyside.png")


# ── Plot 4: Per-sample classified read gain ───────────────────────────────────

def plot_classified_gain(pre, post, out_dir, label="PanHuman"):
    """Show how many more reads are classified post-depletion per sample."""
    df = pd.DataFrame({
        "sample_id":    pre.index,
        "pre_class":    pre["classified_reads"].values,
        "post_class":   post["classified_reads"].values,
        "pre_human":    pre["human_reads"].values,
        "ancestry":     pre["ancestry"].fillna("Unknown").values,
        "total_pre":    pre["total_reads"].values,
    })
    df["nonhuman_pre"]  = df["pre_class"]  - df["pre_human"]
    df["nonhuman_post"] = df["post_class"]
    df["gain_abs"]      = df["nonhuman_post"] - df["nonhuman_pre"]
    df["gain_pct"]      = df["gain_abs"] / df["total_pre"].replace(0, np.nan) * 100
    df = df.sort_values("gain_pct", ascending=False).reset_index(drop=True)

    fig, axes = plt.subplots(1, 2, figsize=(16, 5))
    fig.suptitle(f"Non-Human Classified Read Gain After {label} Depletion",
                 fontsize=13, fontweight="bold")

    x = np.arange(len(df))
    colors = [ANCESTRY_COLORS.get(a, "#757575") for a in df["ancestry"]]

    # Absolute gain
    axes[0].bar(x, df["gain_abs"] / 1e3, color=colors, width=1.0)
    axes[0].axhline(0, color="black", lw=0.8)
    axes[0].set_ylabel("Gain in Non-Human Classified Reads (thousands)")
    axes[0].set_title("Absolute Gain (thousands of reads)")
    axes[0].set_xticks([])
    axes[0].set_xlabel(f"Samples (n={len(df)}, sorted by % gain)")

    # % gain
    axes[1].bar(x, df["gain_pct"], color=colors, width=1.0)
    axes[1].axhline(0, color="black", lw=0.8)
    axes[1].set_ylabel("Gain in Non-Human Classified Reads (% of total)")
    axes[1].set_title("Relative Gain (% of total reads)")
    axes[1].set_xticks([])
    axes[1].set_xlabel(f"Samples (n={len(df)}, sorted by % gain)")
    axes[1].legend(handles=ancestry_legend(), loc="upper right", fontsize=8,
                   title="Ancestry")

    plt.tight_layout()
    plt.savefig(out_dir / "C04_classified_gain.png", dpi=150, bbox_inches="tight")
    plt.close(); print("[PLOT] C04_classified_gain.png")


# ── Plot 5: Shannon delta distribution by ancestry ────────────────────────────

def plot_delta_shannon(pre, post, out_dir, label="PanHuman"):
    anc_present = [a for a in ANC_ORDER if a in pre["ancestry"].values]

    delta_shannon = post["shannon"] - pre["shannon"]
    delta_genera  = post["observed_genera"] - pre["observed_genera"]

    df = pd.DataFrame({
        "sample_id":    pre.index,
        "delta_shannon": delta_shannon.values,
        "delta_genera":  delta_genera.values,
        "ancestry":      pre["ancestry"].fillna("Unknown").values,
    })

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    fig.suptitle(f"Change in Alpha Diversity After {label} Depletion (Post − Pre)",
                 fontsize=13, fontweight="bold")

    for ax, col, label in zip(
        axes,
        ["delta_shannon", "delta_genera"],
        ["ΔShannon Entropy", "ΔObserved Genera"],
    ):
        data_by_anc = [df.loc[df["ancestry"] == a, col].values for a in anc_present]
        bp = ax.boxplot(data_by_anc, patch_artist=True,
                        medianprops=dict(color="black", lw=2))
        for patch, anc in zip(bp["boxes"], anc_present):
            patch.set_facecolor(ANCESTRY_COLORS.get(anc, "#757575"))
            patch.set_alpha(0.8)

        for i, (anc, vals) in enumerate(zip(anc_present, data_by_anc)):
            jitter = np.random.default_rng(42).uniform(-0.15, 0.15, len(vals))
            ax.scatter(i + 1 + jitter, vals,
                       color=ANCESTRY_COLORS.get(anc, "#757575"),
                       alpha=0.45, s=15, zorder=3)

        ax.axhline(0, color="black", ls="--", lw=1, alpha=0.6)
        ax.set_xticklabels(anc_present, rotation=20, ha="right", fontsize=9)
        ax.set_ylabel(label)
        ax.set_title(label)

        # Annotate medians
        for i, vals in enumerate(data_by_anc):
            if len(vals) > 0:
                ax.text(i + 1, np.median(vals), f" {np.median(vals):+.2f}",
                        va="center", fontsize=7, color="black")

    # ── Pairwise post-hoc annotations on ΔShannon panel ──────────────────────
    ax_shannon = axes[0]
    posthoc = pairwise_delta_posthoc(df, "delta_shannon", ANC_ORDER)
    if not posthoc.empty:
        sig_pairs = posthoc[posthoc["significant"]]
        # Only annotate pairs involving ≥10 samples on each side
        sig_pairs = sig_pairs[(sig_pairs["n1"] >= 10) & (sig_pairs["n2"] >= 10)]
        y_max = df["delta_shannon"].max()
        y_step = (y_max - df["delta_shannon"].min()) * 0.10

        for row_i, (_, row) in enumerate(sig_pairs.iterrows()):
            g1_idx = anc_present.index(row["group1"]) + 1
            g2_idx = anc_present.index(row["group2"]) + 1
            y = y_max + y_step * (row_i + 1)
            ax_shannon.plot([g1_idx, g1_idx, g2_idx, g2_idx],
                            [y - y_step*0.3, y, y, y - y_step*0.3],
                            color="black", lw=1.0)
            p_label = (f"p={row['p_adj_BH']:.2e}" if row['p_adj_BH'] >= 0.001
                       else f"p={row['p_adj_BH']:.2e}")
            ax_shannon.text((g1_idx + g2_idx) / 2, y + y_step * 0.05,
                            p_label, ha="center", va="bottom", fontsize=7.5,
                            fontweight="bold")
        # Expand y-axis to fit brackets
        if len(sig_pairs):
            ax_shannon.set_ylim(top=y_max + y_step * (len(sig_pairs) + 2))

        ax_shannon.set_title("\u0394Shannon Entropy\n(brackets = significant pairwise, BH-corrected)",
                             fontsize=9)

    plt.tight_layout()
    plt.savefig(out_dir / "C05_delta_diversity.png", dpi=150, bbox_inches="tight")
    plt.close(); print("[PLOT] C05_delta_diversity.png")


# ── Plot 6: Fixed flagged heatmap with proper per-sample colour bar ───────────

def plot_fixed_heatmap(result_dir, label, out_dir, out_name, cov_threshold=0.005):
    """Re-draw the flagged genera heatmap with a correct per-sample ancestry strip."""
    d = Path(result_dir)
    cov_f = d / "flagged_cov_pivot.csv"
    summ_f = d / "per_sample_summary.csv"
    if not cov_f.exists():
        print(f"[SKIP] {out_name} — flagged_cov_pivot.csv not found"); return

    pivot = pd.read_csv(cov_f, index_col=0)
    summ  = pd.read_csv(summ_f)[["sample_id","ancestry"]].drop_duplicates()
    anc_map = summ.set_index("sample_id")["ancestry"].fillna("Unknown")

    # Sort columns by prevalence
    col_order = (pivot > 0).sum().sort_values(ascending=False).index
    pivot = pivot[col_order]

    # Sort rows by ancestry
    pivot["_anc"] = pivot.index.map(anc_map).fillna("Unknown")
    pivot = pivot.sort_values("_anc").drop(columns="_anc")

    n_rows = len(pivot)
    n_cols = len(pivot.columns)
    fig_h = max(7, n_rows * 0.12)
    fig_w = max(12, n_cols * 0.45)

    fig = plt.figure(figsize=(fig_w + 1.2, fig_h))
    # GridSpec: [ancestry strip | heatmap]
    gs = gridspec.GridSpec(1, 2, width_ratios=[0.025, 1], wspace=0.01)
    ax_strip = fig.add_subplot(gs[0])
    ax_heat  = fig.add_subplot(gs[1])

    # Draw ancestry colour strip
    for i, sid in enumerate(pivot.index):
        anc   = anc_map.get(sid, "Unknown")
        color = ANCESTRY_COLORS.get(anc, "#757575")
        ax_strip.add_patch(plt.Rectangle((0, i), 1, 1, color=color))
    ax_strip.set_xlim(0, 1); ax_strip.set_ylim(0, n_rows)
    ax_strip.set_xticks([]); ax_strip.set_yticks([])
    ax_strip.set_ylabel("Samples (sorted by ancestry)", fontsize=8)

    # Heatmap
    sns.heatmap(pivot, cmap="YlOrRd", linewidths=0.2, linecolor="#f0f0f0",
                ax=ax_heat, cbar_kws={"label": "Coverage", "shrink": 0.5},
                yticklabels=False, xticklabels=True)
    ax_heat.set_title(f"{label} — High-Confidence Genera (coverage ≥ {cov_threshold})",
                      fontsize=10, fontweight="bold")
    ax_heat.set_xlabel("Genus", fontsize=9)
    ax_heat.set_ylabel("")
    plt.setp(ax_heat.get_xticklabels(), rotation=45, ha="right", fontsize=7)

    # Ancestry legend
    legend_handles = [mpatches.Patch(color=c, label=a)
                      for a, c in ANCESTRY_COLORS.items() if a != "Unknown"]
    ax_heat.legend(handles=legend_handles, loc="upper right", fontsize=7,
                   title="Ancestry", bbox_to_anchor=(1.22, 1.0))

    plt.savefig(out_dir / out_name, dpi=150, bbox_inches="tight")
    plt.close(); print(f"[PLOT] {out_name}")


# ── Plot 7: Summary stats table ───────────────────────────────────────────────

def write_comparison_report(pre, post, out_dir, label="PanHuman"):
    lines = []
    sep = "=" * 68
    lines += [sep, f"HEROIC1K — Pre vs Post-Depletion Comparison Report ({label})", sep, ""]

    anc_present = [a for a in ANC_ORDER if a in pre["ancestry"].values]

    lines.append("── Overall Read-Level Changes " + "─" * 30)
    for col, label in [("human_pct","Human %"), ("unclassified_pct","Unclassified %"),
                        ("classified_pct","Classified %")]:
        pre_m  = pre[col].median()
        post_m = post[col].median()
        p      = wilcoxon_pval(pre[col].values, post[col].values)
        lines.append(f"  {label:<22}  pre={pre_m:7.3f}  post={post_m:7.3f}  "
                     f"Δ={post_m-pre_m:+7.3f}  Wilcoxon p={p:.3e}")
    lines.append("")

    lines.append("── Alpha Diversity Changes " + "─" * 34)
    for col, label in [("shannon","Shannon"), ("observed_genera","Observed genera")]:
        pre_m  = pre[col].median()
        post_m = post[col].median()
        p      = wilcoxon_pval(pre[col].values, post[col].values)
        lines.append(f"  {label:<22}  pre={pre_m:7.3f}  post={post_m:7.3f}  "
                     f"Δ={post_m-pre_m:+7.3f}  Wilcoxon p={p:.3e}")
    lines.append("")

    lines.append("── Shannon Change by Ancestry " + "─" * 31)
    lines.append(f"  {'Ancestry':<22}  {'Pre med':>8}  {'Post med':>9}  "
                 f"{'Delta':>7}  {'MW p':>12}  n")
    lines.append("  " + "-" * 68)
    for anc in anc_present:
        pm = pre.loc[pre["ancestry"]==anc, "ancestry"]
        pre_s  = pre.loc[pre["ancestry"]==anc,  "shannon"].values
        post_s = post.loc[post["ancestry"]==anc, "shannon"].values
        if len(pre_s) < 2: continue
        p = mw_pval(pre_s, post_s)
        lines.append(f"  {anc:<22}  {np.median(pre_s):>8.3f}  {np.median(post_s):>9.3f}  "
                     f"{np.median(post_s)-np.median(pre_s):>+7.3f}  {p:>12.3e}  {len(pre_s)}")
    lines.append("")

    lines.append("── Non-Human Classified Read Gain " + "─" * 26)
    pre_nh  = (pre["classified_reads"]  - pre["human_reads"]).values
    post_nh = post["classified_reads"].values
    gain    = post_nh - pre_nh
    gain_pct = gain / pre["total_reads"].values * 100
    lines.append(f"  Median absolute gain:     {np.median(gain):>12,.0f} reads")
    lines.append(f"  Median % of total reads:  {np.median(gain_pct):>12.3f}%")
    lines.append(f"  Samples with gain > 0:    {(gain > 0).sum()} / {len(gain)}")
    lines.append("")

    # ── Pairwise post-hoc: ΔShannon between ancestry groups ──────────────────
    delta_df = pd.DataFrame({
        "sample_id":     pre.index,
        "delta_shannon": (post["shannon"] - pre["shannon"]).values,
        "delta_genera":  (post["observed_genera"] - pre["observed_genera"]).values,
        "ancestry":      pre["ancestry"].fillna("Unknown").values,
    })

    lines.append("── Pairwise ΔShannon Post-Hoc Tests (Mann-Whitney U, BH-corrected) " + "─" * 2)
    lines.append("   This answers: does the AMOUNT of diversity gain differ by ancestry?")
    lines.append("")
    ph = pairwise_delta_posthoc(delta_df, "delta_shannon", ANC_ORDER)
    if not ph.empty:
        lines.append(f"  {'Comparison':<38}  {'Med Δ1':>7}  {'Med Δ2':>7}  "
                     f"{'p_raw':>10}  {'p_adj(BH)':>10}  Sig?")
        lines.append("  " + "-" * 82)
        for _, row in ph.iterrows():
            comp = f"{row['group1']} (n={row['n1']}) vs {row['group2']} (n={row['n2']})"
            sig  = "***" if row["p_adj_BH"] < 0.001 else "**" if row["p_adj_BH"] < 0.01                    else "*" if row["p_adj_BH"] < 0.05 else "ns"
            lines.append(f"  {comp:<38}  {row['median1']:>+7.3f}  {row['median2']:>+7.3f}  "
                         f"{row['p_raw']:>10.3e}  {row['p_adj_BH']:>10.3e}  {sig}")
        lines.append("")

        # Save as CSV too
        ph.to_csv(out_dir / "posthoc_delta_shannon.csv", index=False)
        print("[OUT] posthoc_delta_shannon.csv")

        # Pull out the headline comparison for the paper
        sa_eu = ph[(ph["group1"].str.contains("Southern")) &
                   (ph["group2"].str.contains("European")) |
                   (ph["group1"].str.contains("European")) &
                   (ph["group2"].str.contains("Southern"))]
        if not sa_eu.empty:
            r = sa_eu.iloc[0]
            lines.append("  ── HEADLINE STATISTIC FOR PAPER ──────────────────────────────────")
            lines.append(f"  Southern African ΔShannon median: {r['median1']:+.3f}")
            lines.append(f"  European         ΔShannon median: {r['median2']:+.3f}")
            lines.append(f"  Mann-Whitney U p (raw):           {r['p_raw']:.3e}")
            lines.append(f"  Mann-Whitney U p (BH-adjusted):   {r['p_adj_BH']:.3e}")
            lines.append(f"  Interpretation: European samples gained {abs(r['median2']-r['median1']):.3f}")
            lines.append(f"  more Shannon units post-depletion than Southern African samples,")
            lines.append(f"  consistent with systematic underrepresentation of African genomic")
            lines.append(f"  diversity in {label}-only reference databases.")
            lines.append("")
    else:
        lines.append("  (insufficient group sizes for pairwise tests)")
        lines.append("")

    lines += ["", sep]
    text = "\n".join(lines)
    print(text)
    with open(out_dir / "comparison_report.txt", "w") as fh:
        fh.write(text + "\n")
    print("\n[OUT] comparison_report.txt")


# ── Main ───────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--no_dep_dir", required=True,
                    help="Directory with no-depletion result CSVs")
    ap.add_argument("--dep_dir",    required=True,
                    help="Directory with post-depletion result CSVs")
    ap.add_argument("--output_dir", default="./heroic1k_comparison")
    ap.add_argument("--label",      default="T2T",
                    help="Depletion method label used in plot titles (e.g. T2T, GRCh38, PanHuman)")
    args = ap.parse_args()

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    print(f"[INFO] Output: {out_dir.resolve()}")

    pre_raw  = load_dir(args.no_dep_dir)
    post_raw = load_dir(args.dep_dir)
    pre, post = align(pre_raw, post_raw)

    # Comparison plots
    plot_human_waterfall(pre, post, out_dir, args.label)
    plot_shannon_paired(pre, post, out_dir, args.label)
    plot_ancestry_sidebyside(pre, post, out_dir, args.label)
    plot_classified_gain(pre, post, out_dir, args.label)
    plot_delta_shannon(pre, post, out_dir, args.label)

    # Fixed heatmaps (corrected colour bar)
    plot_fixed_heatmap(args.no_dep_dir, "Pre-Depletion",
                       out_dir, "C06_heatmap_pre_fixed.png")
    plot_fixed_heatmap(args.dep_dir,    "Post-Depletion",
                       out_dir, "C07_heatmap_post_fixed.png")

    write_comparison_report(pre, post, out_dir, args.label)

    print(f"\n[DONE] Results in: {out_dir.resolve()}")


if __name__ == "__main__":
    main()
