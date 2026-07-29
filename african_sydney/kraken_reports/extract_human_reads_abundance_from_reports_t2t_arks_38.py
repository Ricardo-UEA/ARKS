#!/usr/bin/env python3
"""
Extract human read statistics from KrakenUniq reports for three depletion conditions.
Run from: ~/scratch/gen_kmers/data/african_sydney/kraken_reports/

Usage:
    python3 extract_human_reads.py

Output:
    human_reads_summary.csv  — one row per sample per condition
"""

import os
import re
import csv

# ── Config ────────────────────────────────────────────────────────────────────
REPORT_DIRS = {
    "No_depletion":  "no_depletion_kraken_reports",          # fill in if you have a no-depletion report dir
    "GRCh38":        "GRCh38_mcf_05_depletion_reports",
    "T2T":           "t2t_mcf_05_depletion_reports",
    "ARKS":          "ARKS_mcf_05_depletion_reports",
}

# Human taxon identifiers to look for
HUMAN_TAXIDS  = {"9606"}          # Homo sapiens
HUMAN_NAMES   = {"homo sapiens"}  # fallback name match

OUTPUT_FILE = "human_reads_summary.csv"

# ── Parser ────────────────────────────────────────────────────────────────────
def parse_report(filepath):
    """
    Parse a KrakenUniq report and return:
      total_reads      — total reads in sample (classified + unclassified)
      human_reads      — reads assigned to Homo sapiens (taxReads at species level)
      human_pct        — % of total reads
      unclassified_pct — % unclassified
    """
    total_reads       = None
    human_reads       = 0
    human_pct         = 0.0
    unclassified_reads = 0

    try:
        with open(filepath) as f:
            for line in f:
                if line.startswith("#"):
                    continue

                parts = line.rstrip("\n").split("\t")
                if len(parts) < 9:
                    continue

                pct_str   = parts[0].strip()
                reads_str = parts[1].strip()
                taxreads  = parts[2].strip()
                rank      = parts[7].strip()
                taxname   = parts[8].strip().lower()
                taxid     = parts[6].strip()

                try:
                    reads = int(reads_str)
                    pct   = float(pct_str)
                except ValueError:
                    continue

                # Total reads = unclassified + root (root covers everything classified)
                # Easiest: unclassified row has taxID 0
                if taxid == "0":
                    unclassified_reads = reads
                    continue

                # Root row (taxID 1) gives total classified
                if taxid == "1" and total_reads is None:
                    # total = unclassified + classified
                    # We'll compute after
                    classified_reads = reads

                # Human: Homo sapiens species row
                if taxid in HUMAN_TAXIDS or "homo sapiens" in taxname:
                    if rank in ("species", "no rank") or "homo sapiens" == taxname:
                        try:
                            human_reads = int(taxreads)  # taxReads = reads at this node only
                            human_pct   = pct
                        except ValueError:
                            pass

        # Compute total
        if unclassified_reads > 0:
            total_reads = unclassified_reads + classified_reads

    except Exception as e:
        print(f"  ERROR parsing {filepath}: {e}")
        return None

    if total_reads is None or total_reads == 0:
        return None

    # Recompute human % from raw counts for accuracy
    human_pct_computed = (human_reads / total_reads) * 100.0

    return {
        "total_reads":    total_reads,
        "human_reads":    human_reads,
        "human_pct":      human_pct_computed,
        "unclassified":   unclassified_reads,
    }


# ── Main ──────────────────────────────────────────────────────────────────────
def main():
    rows = []

    for condition, report_dir in REPORT_DIRS.items():
        if report_dir is None:
            print(f"Skipping {condition} — no directory set")
            continue

        if not os.path.isdir(report_dir):
            print(f"WARNING: {report_dir} not found — skipping {condition}")
            continue

        report_files = sorted([
            f for f in os.listdir(report_dir)
            if f.endswith(".report") or f.endswith(".txt") or "krakenuniq" in f.lower()
        ])

        print(f"{condition}: found {len(report_files)} report files in {report_dir}")

        for fname in report_files:
            # Extract sample ID from filename (everything before first dot)
            sample_id = fname.split(".")[0]
            fpath = os.path.join(report_dir, fname)
            result = parse_report(fpath)

            if result is None:
                print(f"  Could not parse {fname}")
                continue

            rows.append({
                "sample_id":    sample_id,
                "condition":    condition,
                "total_reads":  result["total_reads"],
                "human_reads":  result["human_reads"],
                "human_pct":    round(result["human_pct"], 6),
                "unclassified": result["unclassified"],
            })

    # Write output
    if rows:
        fieldnames = ["sample_id", "condition", "total_reads",
                      "human_reads", "human_pct", "unclassified"]
        with open(OUTPUT_FILE, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(rows)
        print(f"\nDone. {len(rows)} rows written to {OUTPUT_FILE}")

        # Quick summary
        from collections import defaultdict
        import statistics
        by_cond = defaultdict(list)
        for r in rows:
            by_cond[r["condition"]].append(r["human_pct"])
        print("\nQuick summary (human read %):")
        for cond, vals in by_cond.items():
            print(f"  {cond:15s}  n={len(vals):3d}  "
                  f"median={statistics.median(vals):.4f}%  "
                  f"mean={statistics.mean(vals):.4f}%")
    else:
        print("No data extracted — check report directory paths")


if __name__ == "__main__":
    main()
