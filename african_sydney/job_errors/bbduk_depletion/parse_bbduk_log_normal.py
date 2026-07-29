import re
import sys
import csv
import os


def parse_bbduk_log(log_file):
    with open(log_file, 'r', encoding='utf-8', errors='ignore') as f:
        log = f.read()

    result = {
        'sample': 'NA',
        'reference': 'NA',
        'input_reads': 'NA',
        'input_bases': 'NA',
        'removed_reads': 'NA',
        'removed_bases': 'NA',
        'remaining_reads': 'NA',
        'remaining_bases': 'NA',
        'read_speed': 'NA',
        'base_speed': 'NA'
    }

    # -----------------------------
    # Sample name
    # Try Snakemake wildcards first
    # -----------------------------
    sample_match = re.search(r"wildcards:\s+sample=([^\s]+)", log)
    if sample_match:
        result['sample'] = sample_match.group(1)
    else:
        # Try to extract from in1=.../10651_R1.fastq.gz
        in1_match = re.search(r'in1=([^\s,]+)', log)
        if in1_match:
            in1_path = in1_match.group(1)
            filename = os.path.basename(in1_path)

            # Remove common paired-end suffixes
            sample = re.sub(r'(_R?1(?:_001)?|\.R1)\.f(ast)?q(\.gz)?$', '', filename, flags=re.IGNORECASE)
            result['sample'] = sample
        else:
            # Fallback: use the log filename
            result['sample'] = os.path.basename(log_file).split('.')[0]

    # -----------------------------
    # Reference
    # -----------------------------
    ref_match = re.search(r"ref=([^\s,]+)", log)
    if ref_match:
        result['reference'] = os.path.basename(ref_match.group(1))

    # -----------------------------
    # Input reads and bases
    # Example:
    # Input: 2256114 reads 333345223 bases.
    # -----------------------------
    input_match = re.search(
        r"Input:\s+([\d,]+)\s+reads\s+([\d,]+)\s+bases",
        log
    )
    if input_match:
        result['input_reads'] = int(input_match.group(1).replace(",", ""))
        result['input_bases'] = int(input_match.group(2).replace(",", ""))

    # -----------------------------
    # Removed reads and bases
    # Example:
    # Total Removed: 73754 reads (3.27%) 10570055 bases (3.17%)
    # -----------------------------
    removed_match = re.search(
        r"Total Removed:\s+([\d,]+)\s+reads.*?\s+([\d,]+)\s+bases",
        log,
        flags=re.DOTALL
    )
    if removed_match:
        result['removed_reads'] = int(removed_match.group(1).replace(",", ""))
        result['removed_bases'] = int(removed_match.group(2).replace(",", ""))

    # -----------------------------
    # Remaining reads and bases
    # Example:
    # Result: 2182360 reads (96.73%) 322775168 bases (96.83%)
    # -----------------------------
    remaining_match = re.search(
        r"Result:\s+([\d,]+)\s+reads.*?\s+([\d,]+)\s+bases",
        log,
        flags=re.DOTALL
    )
    if remaining_match:
        result['remaining_reads'] = int(remaining_match.group(1).replace(",", ""))
        result['remaining_bases'] = int(remaining_match.group(2).replace(",", ""))

    # -----------------------------
    # Read speed
    # Example:
    # Reads Processed: 2256k 5.58k reads/sec
    # -----------------------------
    read_speed_match = re.search(
        r"Reads Processed:\s+\S+\s+([\d.]+[kKmMgG]?)\s+reads/sec",
        log
    )
    if read_speed_match:
        result['read_speed'] = read_speed_match.group(1)

    # -----------------------------
    # Base speed
    # Example:
    # Bases Processed: 333m 0.82m bases/sec
    # -----------------------------
    base_speed_match = re.search(
        r"Bases Processed:\s+\S+\s+([\d.]+[kKmMgG]?)\s+bases/sec",
        log
    )
    if base_speed_match:
        result['base_speed'] = base_speed_match.group(1)

    return result


def main(log_files):
    fieldnames = [
        'sample', 'reference', 'input_reads', 'input_bases',
        'removed_reads', 'removed_bases',
        'remaining_reads', 'remaining_bases',
        'read_speed', 'base_speed'
    ]

    writer = csv.DictWriter(sys.stdout, fieldnames=fieldnames)
    writer.writeheader()

    for log_file in log_files:
        result = parse_bbduk_log(log_file)
        writer.writerow(result)


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python parse_bbduk_log.py <log_file1> <log_file2> ...", file=sys.stderr)
        sys.exit(1)

    main(sys.argv[1:])
