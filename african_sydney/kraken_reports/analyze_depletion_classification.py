#!/usr/bin/env python3
"""
analyze_depletion_classification.py

Compares KrakenUniq classification reports before vs after host depletion
(no-depletion baseline vs ARKS MCF=0.5 depleted reads) to identify which
microbial taxa are:

  (a) intrinsically hard for KrakenUniq to classify with confidence
      (low unique k-mer coverage / high k-mer duplication against the DB,
      regardless of depletion),
  (b) disproportionately lost after depletion (candidate off-target
      removal by the host k-mer set), or
  (c) robustly and confidently detected in both conditions.

Input: two directories of KrakenUniq --report-file outputs, one file per
sample, matched by filename (or by a normalised sample ID if filenames
differ between the two directories -- see --strip-patterns).

KrakenUniq report columns (tab-separated; a header row is auto-detected
and skipped if present):
  1. pct        % of reads in the clade rooted at this taxon
  2. reads      reads in the clade (this taxon + descendants)
  3. taxReads   reads assigned directly to this taxon
  4. kmers      unique k-mers assigned to this taxon/clade
  5. dup        avg. number of times each unique k-mer was seen
                (duplication -- higher means more repetitive/ambiguous)
  6. cov        coverage of this taxon's unique k-mers in the DB
                (closer to 1 = confident; closer to 0 = poorly covered,
                often driven by shared k-mers with relatives or gaps in
                the reference)
  7. taxID
  8. rank       full rank name as written by KrakenUniq (e.g. "no rank",
                "superkingdom", "phylum", "family", "genus", "species",
                "subspecies", "strain", "species group", ...) -- NOT the
                single-letter codes used by plain Kraken/Kraken2 reports
  9. taxName    scientific name (indented by tree depth with spaces)

Real KrakenUniq report files also start with 1-2 "#"-prefixed metadata
lines (version, DB path, full command line) and then a tab-separated
header row ("%	reads	taxReads	kmers	dup	cov	taxID	rank	taxName")
before the data rows -- all of which this script strips automatically.

Usage:
  python3 analyze_depletion_classification.py \
      --before /gpfs/home/hce24xau/scratch/gen_kmers/data/african_sydney/kraken_reports/no_depletion_kraken_reports \
      --after  /gpfs/home/hce24xau/scratch/gen_kmers/data/african_sydney/kraken_reports/arks_mcf_05_depletion_reports \
      --rank genus \
      --out-prefix heroic1k_african_sydney

Outputs (written to --out-dir, default: current directory):
  {prefix}_long.csv            per-sample per-taxon rows, both conditions
  {prefix}_taxon_summary.csv   per-taxon aggregate stats + before/after deltas + flags
  stdout                       quick top-N summaries + sample-matching diagnostics

Notes:
  - This is a first-pass heuristic screen, not a validated statistical
    test. The --cov-low / --dup-high / --detection-drop / --reads-log2fc-drop
    thresholds are deliberately exposed as CLI args so they can be tuned
    and the sensitivity of the categorisation checked before anything is
    written up.
  - No third-party dependencies (stdlib only), so it should run on Ada
    without needing a particular Python environment/module loaded.
"""

import argparse
import csv
import math
import os
import re
import statistics
import sys
from collections import defaultdict

DEFAULT_STRIP_PATTERNS = [
    r"\.kreport$", r"\.report$", r"\.txt$", r"\.tsv$", r"\.gz$",
    r"_no_?depletion", r"_raw", r"_baseline",
    r"_arks(_mcf_?0?5)?", r"_depleted", r"_depletion",
]


def extract_sample_id(filename, strip_patterns):
    sid = filename
    for pat in strip_patterns:
        sid = re.sub(pat, "", sid, flags=re.IGNORECASE)
    return sid.strip("._-")


def parse_report(path):
    """Yield dicts, one per taxon row, for a single KrakenUniq report file."""
    with open(path, "r", errors="replace") as fh:
        lines = [ln.rstrip("\n") for ln in fh if ln.strip()]
    # Drop KrakenUniq's "#"-prefixed metadata lines (version, DB, command line)
    lines = [ln for ln in lines if not ln.startswith("#")]
    if not lines:
        return
    # Auto-detect + skip the tab-separated header row if present (first field isn't numeric)
    start = 0
    first_fields = lines[0].split("\t")
    try:
        float(first_fields[0])
    except (ValueError, IndexError):
        start = 1
    for ln in lines[start:]:
        fields = ln.split("\t")
        if len(fields) < 9:
            continue
        pct, reads, taxreads, kmers, dup, cov, taxid, rank, taxname = fields[:9]
        try:
            yield {
                "pct": float(pct),
                "reads": int(reads),
                "taxReads": int(taxreads),
                "kmers": int(kmers),
                "dup": float(dup),
                "cov": float(cov),
                "taxID": taxid.strip(),
                "rank": rank.strip(),
                "taxName": taxname.strip(),
            }
        except ValueError:
            # malformed row -- skip rather than crash the whole run
            continue


def collect_samples(directory):
    files = {}
    for fn in os.listdir(directory):
        full = os.path.join(directory, fn)
        if os.path.isfile(full):
            files[fn] = full
    return files


def match_samples(before_files, after_files, strip_patterns):
    """Return dict sample_id -> (before_path, after_path), plus unmatched lists."""
    before_by_id = {extract_sample_id(fn, strip_patterns): p for fn, p in before_files.items()}
    after_by_id = {extract_sample_id(fn, strip_patterns): p for fn, p in after_files.items()}

    matched = {}
    for sid, bpath in before_by_id.items():
        if sid in after_by_id:
            matched[sid] = (bpath, after_by_id[sid])

    unmatched_before = sorted(set(before_by_id) - set(after_by_id))
    unmatched_after = sorted(set(after_by_id) - set(before_by_id))
    return matched, unmatched_before, unmatched_after


def median(vals):
    return statistics.median(vals) if vals else float("nan")


def log2fc(after, before, pseudocount=1.0):
    return math.log2((after + pseudocount) / (before + pseudocount))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--before", default="/gpfs/home/hce24xau/scratch/gen_kmers/data/african_sydney/kraken_reports/no_depletion_kraken_reports",
                     help="Directory of no-depletion KrakenUniq reports")
    ap.add_argument("--after", default="/gpfs/home/hce24xau/scratch/gen_kmers/data/african_sydney/kraken_reports/arks_mcf_05_depletion_reports",
                     help="Directory of post-depletion (e.g. ARKS MCF=0.5) KrakenUniq reports")
    ap.add_argument("--rank", default="genus", choices=["genus", "species", "G", "S", "both"],
                     help="Taxonomic rank to analyse -- matches KrakenUniq's full rank-name strings, "
                          "not single-letter codes (default: genus)")
    ap.add_argument("--min-reads", type=int, default=10, help="Minimum clade reads to count a taxon as 'detected' in a sample (default: 10)")
    ap.add_argument("--cov-low", type=float, default=0.1, help="cov below this = poor k-mer coverage / low confidence (default: 0.1)")
    ap.add_argument("--dup-high", type=float, default=5.0, help="dup above this = high k-mer duplication / ambiguous (default: 5.0)")
    ap.add_argument("--detection-drop", type=float, default=0.2, help="Drop in detection rate (before->after) to flag as depletion-sensitive (default: 0.2 = 20 percentage points)")
    ap.add_argument("--reads-log2fc-drop", type=float, default=-1.0, help="log2 fold-change in median reads (after vs before) below this = flag as depletion-sensitive (default: -1.0 = >=50%% reduction)")
    ap.add_argument("--strip-patterns", default=None, help="Comma-separated extra regex patterns to strip from filenames when matching samples between directories")
    ap.add_argument("--out-prefix", default="classification_comparison", help="Prefix for output CSV files")
    ap.add_argument("--out-dir", default=".", help="Directory to write output CSVs")
    args = ap.parse_args()

    strip_patterns = list(DEFAULT_STRIP_PATTERNS)
    if args.strip_patterns:
        strip_patterns = [p.strip() for p in args.strip_patterns.split(",")] + strip_patterns

    rank_alias = {"g": "genus", "genus": "genus", "s": "species", "species": "species"}
    if args.rank == "both":
        ranks = {"genus", "species"}
    else:
        ranks = {rank_alias[args.rank.lower()]}

    before_files = collect_samples(args.before)
    after_files = collect_samples(args.after)
    print(f"[info] {len(before_files)} files in --before, {len(after_files)} files in --after", file=sys.stderr)

    matched, unmatched_before, unmatched_after = match_samples(before_files, after_files, strip_patterns)
    print(f"[info] matched {len(matched)} samples across both conditions", file=sys.stderr)
    if unmatched_before:
        print(f"[warn] {len(unmatched_before)} 'before' samples had no match in 'after' (first 5: {unmatched_before[:5]})", file=sys.stderr)
    if unmatched_after:
        print(f"[warn] {len(unmatched_after)} 'after' samples had no match in 'before' (first 5: {unmatched_after[:5]})", file=sys.stderr)
    if not matched:
        print("[error] no samples matched between --before and --after directories.", file=sys.stderr)
        print("        Inspect filenames below and pass --strip-patterns to normalise them.", file=sys.stderr)
        print("        before filenames (sample):", list(before_files)[:10], file=sys.stderr)
        print("        after filenames (sample):", list(after_files)[:10], file=sys.stderr)
        sys.exit(1)

    os.makedirs(args.out_dir, exist_ok=True)
    long_path = os.path.join(args.out_dir, f"{args.out_prefix}_long.csv")
    summary_path = os.path.join(args.out_dir, f"{args.out_prefix}_taxon_summary.csv")

    # per-taxon accumulator: (taxID, rank) -> condition -> list of per-sample metric rows
    per_taxon = defaultdict(lambda: {"before": [], "after": []})
    taxon_meta = {}  # (taxID, rank) -> taxName

    with open(long_path, "w", newline="") as long_fh:
        long_writer = csv.writer(long_fh)
        long_writer.writerow(["sample_id", "condition", "taxID", "rank", "taxName", "reads", "taxReads", "kmers", "dup", "cov", "pct"])

        for sid, (bpath, apath) in sorted(matched.items()):
            for condition, path in (("before", bpath), ("after", apath)):
                for row in parse_report(path):
                    if row["rank"].lower() not in ranks:
                        continue
                    long_writer.writerow([sid, condition, row["taxID"], row["rank"], row["taxName"],
                                           row["reads"], row["taxReads"], row["kmers"], row["dup"], row["cov"], row["pct"]])
                    key = (row["taxID"], row["rank"])
                    taxon_meta[key] = row["taxName"]
                    per_taxon[key][condition].append(row)

    n_samples = len(matched)
    summary_rows = []
    for key, by_cond in per_taxon.items():
        rank = key[1]
        taxname = taxon_meta[key]
        b = by_cond["before"]
        a = by_cond["after"]

        b_detected = [r for r in b if r["reads"] >= args.min_reads]
        a_detected = [r for r in a if r["reads"] >= args.min_reads]

        b_det_rate = len(b_detected) / n_samples
        a_det_rate = len(a_detected) / n_samples

        b_med_reads = median([r["reads"] for r in b_detected])
        a_med_reads = median([r["reads"] for r in a_detected])
        b_med_cov = median([r["cov"] for r in b_detected])
        a_med_cov = median([r["cov"] for r in a_detected])
        b_med_dup = median([r["dup"] for r in b_detected])
        a_med_dup = median([r["dup"] for r in a_detected])

        reads_l2fc = (
            log2fc(a_med_reads if a_detected else 0.0, b_med_reads if b_detected else 0.0)
            if (a_detected or b_detected) else float("nan")
        )
        det_rate_change = a_det_rate - b_det_rate

        intrinsically_ambiguous = bool(b_detected) and (
            (not math.isnan(b_med_cov) and b_med_cov < args.cov_low) or
            (not math.isnan(b_med_dup) and b_med_dup > args.dup_high)
        )
        depletion_sensitive = bool(b_detected) and (
            det_rate_change <= -args.detection_drop or
            (not math.isnan(reads_l2fc) and reads_l2fc <= args.reads_log2fc_drop)
        )
        is_robust = (
            bool(b_detected) and bool(a_detected)
            and not intrinsically_ambiguous and not depletion_sensitive
            and det_rate_change > -0.05
        )

        summary_rows.append({
            "taxID": key[0], "rank": rank, "taxName": taxname,
            "n_samples": n_samples,
            "before_detection_rate": round(b_det_rate, 3),
            "after_detection_rate": round(a_det_rate, 3),
            "detection_rate_change": round(det_rate_change, 3),
            "before_median_reads": b_med_reads if not math.isnan(b_med_reads) else "",
            "after_median_reads": a_med_reads if not math.isnan(a_med_reads) else "",
            "reads_log2fc": round(reads_l2fc, 3) if not math.isnan(reads_l2fc) else "",
            "before_median_cov": round(b_med_cov, 4) if not math.isnan(b_med_cov) else "",
            "after_median_cov": round(a_med_cov, 4) if not math.isnan(a_med_cov) else "",
            "before_median_dup": round(b_med_dup, 3) if not math.isnan(b_med_dup) else "",
            "after_median_dup": round(a_med_dup, 3) if not math.isnan(a_med_dup) else "",
            "intrinsically_ambiguous": intrinsically_ambiguous,
            "depletion_sensitive": depletion_sensitive,
            "robust": is_robust,
        })

    summary_rows.sort(key=lambda r: (r["before_detection_rate"], r["taxName"]), reverse=True)
    fieldnames = list(summary_rows[0].keys()) if summary_rows else []
    with open(summary_path, "w", newline="") as sfh:
        writer = csv.DictWriter(sfh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(summary_rows)

    print(f"\n[done] wrote {long_path}")
    print(f"[done] wrote {summary_path}\n")

    def show(rows, label, n=15):
        print(f"--- {label} (top {n}) ---")
        for r in rows[:n]:
            print(f"  {r['taxName']:<35} rank={r['rank']:<2} before_det={r['before_detection_rate']:.2f} "
                  f"after_det={r['after_detection_rate']:.2f} log2FC={r['reads_log2fc']!s:<7} "
                  f"cov(before/after)={r['before_median_cov']}/{r['after_median_cov']} "
                  f"dup(before/after)={r['before_median_dup']}/{r['after_median_dup']}")
        print()

    ambiguous = [r for r in summary_rows if r["intrinsically_ambiguous"]]
    sensitive = [r for r in summary_rows if r["depletion_sensitive"]]
    robust_rows = [r for r in summary_rows if r["robust"]]

    ambiguous.sort(key=lambda r: (r["before_median_cov"] if r["before_median_cov"] != "" else 1.0))
    sensitive.sort(key=lambda r: (r["reads_log2fc"] if r["reads_log2fc"] != "" else 0))
    robust_rows.sort(key=lambda r: -r["before_detection_rate"])

    show(ambiguous, "Intrinsically hard to classify (low k-mer coverage / high duplication, both conditions)")
    show(sensitive, "Depletion-sensitive (disproportionately lost after ARKS depletion)")
    show(robust_rows, "Robust / confidently classified in both conditions")


if __name__ == "__main__":
    main()
