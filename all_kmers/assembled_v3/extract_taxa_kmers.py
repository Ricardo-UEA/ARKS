#!/usr/bin/env python3
"""
extract_taxa_kmers.py
---------------------
Extract k-mer sequences from an ARKS FASTA file for a set of target taxa,
based on KrakenUniq classification output.

Usage examples
--------------
# Taxa IDs inline:
python extract_taxa_kmers.py \
    --report classified_kmers_kraken.tsv \
    --fasta  ARKS_master_kmers_filtered.fasta \
    --outdir taxa_kmers_to_remove \
    --taxa 562 87882 1280

# Taxa IDs from a file (one per line):
python extract_taxa_kmers.py \
    --report classified_kmers_kraken.tsv \
    --fasta  ARKS_master_kmers_filtered.fasta \
    --outdir taxa_kmers_to_remove \
    --taxa-file taxa_ids.txt

# Combined (file + extra inline IDs):
python extract_taxa_kmers.py ... --taxa-file taxa_ids.txt --taxa 999 888
"""

import os
import sys
import argparse
from collections import defaultdict


# ------------------------------------------------------------------
# Default taxa IDs (mirrors the original hardcoded list)
# ------------------------------------------------------------------
DEFAULT_TAXA_IDS = [
    "2121", "9", "3050299", "11856", "562", "86661", "2071627", "336810",
    "111527", "2108470", "47715", "87882", "2099", "210", "28450", "114186",
    "1280", "28901", "909768", "573", "717610", "656088", "470", "515350",
    "485", "56", "1282", "1349409", "2653194", "2946593", "2931930", "136841",
    "28903", "45617", "653685", "287", "2107707", "59201", "2371", "2208",
    "643453", "346", "339", "1349410", "1938374", "2838335", "1126", "357276",
    "2653200", "1405", "208962", "2102", "383372", "120962", "1906", "2209",
    "580165", "712361", "3004094", "492670", "1491", "1396", "42879", "171284",
    "1922217", "1428", "669", "337", "1499987", "629295", "54571", "373994",
    "1502", "40477", "181082", "727", "347", "2565304", "2107708", "2782013",
    "2771012", "2294034", "136845", "1653831", "39441", "213615", "44822",
    "2098", "1155739", "83656", "472834", "1444770", "75105", "87883",
    "136849", "1219", "654", "1345117", "1434100", "47760", "3050256", "1351",
    "343", "2597325", "3050269", "1238", "3018339", "68278", "95486", "435897",
    "644357", "1173022", "241425", "1590", "294671", "487", "584", "1254432",
    "249567", "1647413", "747523", "995085", "136842", "747522", "238834",
    "40324", "136843", "2035", "1597781", "46234", "91844", "1520", "372461",
    "1434111", "2839105", "170861", "354276", "412022", "441156", "32630",
    "1401325", "295358", "2879117", "1316932", "262316", "3050254", "743965",
    "363253", "436113", "1249471", "2365036", "561275", "1234378", "1819678",
    "1971437", "2560781", "1311760", "662596"
]
# ------------------------------------------------------------------
# Argument parsing
# ------------------------------------------------------------------
def parse_args():
    parser = argparse.ArgumentParser(
        description="Extract k-mers from an ARKS FASTA by KrakenUniq taxa ID.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "--report", "-r",
        required=True,
        metavar="TSV",
        help="KrakenUniq classified k-mers TSV (classified_kmers_kraken.tsv).",
    )
    parser.add_argument(
        "--fasta", "-f",
        required=True,
        metavar="FASTA",
        help="ARKS master k-mers FASTA file.",
    )
    parser.add_argument(
        "--outdir", "-o",
        default="taxa_kmers_to_remove",
        metavar="DIR",
        help="Output directory (created if absent). Default: taxa_kmers_to_remove",
    )
    parser.add_argument(
        "--taxa", "-t",
        nargs="*",
        metavar="ID",
        default=[],
        help="One or more taxa IDs (space-separated). Combined with --taxa-file.",
    )
    parser.add_argument(
        "--taxa-file",
        metavar="FILE",
        help="Plain-text file with one taxa ID per line (comments with # ignored).",
    )
    parser.add_argument(
        "--use-defaults",
        action="store_true",
        help="Include the built-in default taxa list (used when no --taxa / --taxa-file given).",
    )
    return parser.parse_args()


# ------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------
def load_taxa_ids(args) -> set:
    """Merge taxa IDs from all sources into a single set."""
    ids = set(args.taxa)

    if args.taxa_file:
        with open(args.taxa_file) as fh:
            for line in fh:
                line = line.strip()
                if line and not line.startswith("#"):
                    ids.add(line)

    # Fall back to defaults when nothing was supplied
    if not ids or args.use_defaults:
        ids.update(DEFAULT_TAXA_IDS)

    return ids


def build_kmer_to_taxa_map(report_path: str, target_taxa: set) -> dict:
    """
    Pass 1 — read the report TSV and build a reverse map:
        kmer_id -> set of matching taxa IDs

    This is O(lines) instead of O(lines × taxa).
    """
    kmer_to_taxa = defaultdict(set)
    skipped = 0

    with open(report_path) as fh:
        for line in fh:
            parts = line.split()          # avoids strip+split overhead
            if len(parts) < 3:
                skipped += 1
                continue
            kmer_id = parts[1]
            taxa_id = parts[2]
            if taxa_id in target_taxa:
                kmer_to_taxa[kmer_id].add(taxa_id)

    if skipped:
        print(f"  [warn] skipped {skipped} malformed lines in report.", file=sys.stderr)

    return kmer_to_taxa


def extract_sequences(fasta_path: str, kmer_to_taxa: dict, target_taxa: set) -> dict:
    """
    Pass 2 — stream the FASTA and bucket sequences by taxa.
    Each k-mer is looked up once (O(1)) via the reverse map.
    """
    results = defaultdict(list)

    with open(fasta_path) as fh:
        header = None
        for line in fh:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith(">"):
                header = line[1:]           # strip leading '>'
            elif header is not None:
                matched_taxa = kmer_to_taxa.get(header)
                if matched_taxa:
                    seq = line
                    for taxa_id in matched_taxa:
                        results[taxa_id].append(seq)
                header = None               # reset; ready for next pair

    return results


def write_outputs(output_dir: str, target_taxa: set, results: dict) -> None:
    """Write one file per taxa ID; report counts."""
    os.makedirs(output_dir, exist_ok=True)
    found = 0

    for taxa_id in sorted(target_taxa, key=lambda x: int(x) if x.isdigit() else x):
        kmers = results.get(taxa_id, [])
        out_path = os.path.join(output_dir, f"{taxa_id}_kmers.txt")
        with open(out_path, "w") as out:
            if kmers:
                found += 1
                for kmer in sorted(set(kmers)):   # deduplicate + sort
                    out.write(kmer + "\n")
            else:
                out.write(f"No matching k-mers found for taxa ID {taxa_id}.\n")
        print(f"  [{len(kmers):>7,} k-mers]  {out_path}")

    print(f"\nDone. {found}/{len(target_taxa)} taxa had matching k-mers.")


# ------------------------------------------------------------------
# Main
# ------------------------------------------------------------------
def main():
    args = parse_args()

    # ---- Resolve taxa IDs ----------------------------------------
    target_taxa = load_taxa_ids(args)
    print(f"Target taxa: {len(target_taxa)} IDs")

    # ---- Pass 1: report → kmer→taxa map --------------------------
    print(f"Reading report:  {args.report}")
    kmer_to_taxa = build_kmer_to_taxa_map(args.report, target_taxa)
    print(f"  {len(kmer_to_taxa):,} k-mer IDs matched at least one target taxon.")

    if not kmer_to_taxa:
        print("No k-mers matched any target taxa — nothing to write.", file=sys.stderr)
        sys.exit(0)

    # ---- Pass 2: FASTA → sequences by taxa -----------------------
    print(f"Scanning FASTA:  {args.fasta}")
    results = extract_sequences(args.fasta, kmer_to_taxa, target_taxa)

    # ---- Write outputs -------------------------------------------
    print(f"Writing outputs to: {args.outdir}/")
    write_outputs(args.outdir, target_taxa, results)


if __name__ == "__main__":
    main()
