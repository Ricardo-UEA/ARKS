#!/usr/bin/env python3
"""
Extract metadata/statistics from BBDuk (BBMap) log files into a single CSV.

A BBDuk log looks like the console output of a `jgi.BBDuk` run, e.g.:

    Executing jgi.BBDuk [-Xmx200g, in1=..., in2=..., out=..., ref=..., k=31, ...]
    Version 39.13
    Set INTERLEAVED to true
    ...
    Input:                          2256114 reads           333345223 bases.
    Contaminants:                   30110 reads (1.33%)     4146700 bases (1.24%)
    Total Removed:                  30110 reads (1.33%)     4146700 bases (1.24%)
    Result:                         2226004 reads (98.67%)  329198523 bases (98.76%)
    Time:                           998.181 seconds.
    Reads Processed:       2256k    2.26k reads/sec
    Bases Processed:        333m    0.33m bases/sec

This script walks a directory, finds every file whose contents look like a
BBDuk log (contains "jgi.BBDuk"), pulls out every piece of metadata it can
find, and writes one row per file to a CSV.

Usage:
    python3 extract_bbduk_metadata.py /path/to/logs
    python3 extract_bbduk_metadata.py /path/to/logs -o summary.csv
    python3 extract_bbduk_metadata.py /path/to/logs --no-recursive
    python3 extract_bbduk_metadata.py /path/to/logs --pattern "*_bbduk_stats.txt"
"""
import argparse
import csv
import os
import re
import sys
from fnmatch import fnmatch

# ---------------------------------------------------------------------------
# Regexes for the various pieces of a BBDuk log
# ---------------------------------------------------------------------------
RE_COMMAND_LIST = re.compile(r"Executing\s+jgi\.BBDuk\s*\[(.*)\]\s*$", re.MULTILINE)
RE_VERSION = re.compile(r"^Version\s+(\S+)", re.MULTILINE)
RE_SET = re.compile(r"^Set\s+(\w+)\s+to\s+(\w+)\s*$", re.MULTILINE)
RE_RESET = re.compile(r"^Reset\s+(\w+)\s+to\s+(\w+)(?:\s+because\s+(.*?))?\.?\s*$", re.MULTILINE)
RE_STARTUP_TIME = re.compile(r"^([\d.]+)\s+seconds\.\s*$", re.MULTILINE)
RE_MEMORY = re.compile(
    r"^Memory:\s*max=(\d+)m,\s*total=(\d+)m,\s*free=(\d+)m,\s*used=(\d+)m\s*$",
    re.MULTILINE,
)
RE_KMERS_ADDED = re.compile(
    r"^Added\s+(\d+)\s+kmers;\s*time:\s*([\d.]+)\s+seconds\.\s*$", re.MULTILINE
)
RE_INPUT_MODE = re.compile(r"^Input is being processed as\s+(.+?)\.?\s*$", re.MULTILINE)
RE_OUTPUT_STREAMS = re.compile(r"^Started output streams:\s*([\d.]+)\s+seconds\.\s*$", re.MULTILINE)
RE_PROCESSING_TIME = re.compile(r"^Processing time:\s*([\d.]+)\s+seconds\.\s*$", re.MULTILINE)
RE_TOTAL_TIME = re.compile(r"^Time:\s*([\d.]+)\s+seconds\.\s*$", re.MULTILINE)
RE_RATE_LINE = re.compile(
    r"^(Reads Processed|Bases Processed):\s*(\S+)\s+(\S+)\s+(reads/sec|bases/sec)\s*$",
    re.MULTILINE,
)
# Generic "Label:   N reads (P%)   M bases (Q%)" lines: Input/Contaminants/
# Total Removed/Result/QTrimmed/etc, whatever happens to be present.
RE_STAT_PCT = re.compile(
    r"^([A-Za-z][A-Za-z0-9 /_-]*?):\s+(\d+)\s+reads\s+\(([\d.]+)%\)\s+(\d+)\s+bases\s+\(([\d.]+)%\)\s*$",
    re.MULTILINE,
)
# "Input:   N reads   M bases." (no percentages)
RE_STAT_PLAIN = re.compile(
    r"^([A-Za-z][A-Za-z0-9 /_-]*?):\s+(\d+)\s+reads\s+(\d+)\s+bases\.?\s*$",
    re.MULTILINE,
)

SIGNATURE = "jgi.BBDuk"

# Columns that describe the *run* (paths, CLI flags, tool version) rather than
# a measured result. Everything else found in a log is treated as a "metric".
CONFIG_PREFIXES = ("param_", "set_", "reset_")
CONFIG_EXTRAS = {"jvm_flags", "bbduk_version"}


def is_config_column(name):
    return name.startswith(CONFIG_PREFIXES) or name in CONFIG_EXTRAS


# Preferred left-to-right column order for each output. Anything discovered
# that isn't listed here (e.g. an unusual extra stat line) is appended at the
# end in first-seen order, so nothing is ever silently dropped.
METRIC_ORDER = [
    "input_reads", "input_bases",
    "contaminants_reads", "contaminants_reads_pct", "contaminants_bases", "contaminants_bases_pct",
    "total_removed_reads", "total_removed_reads_pct", "total_removed_bases", "total_removed_bases_pct",
    "result_reads", "result_reads_pct", "result_bases", "result_bases_pct",
    "kmers_added", "kmer_load_time_seconds",
    "startup_time_seconds", "output_stream_start_seconds", "processing_time_seconds", "total_time_seconds",
    "reads_processed_count", "reads_processed_rate", "bases_processed_count", "bases_processed_rate",
    "initial_mem_max_mb", "initial_mem_total_mb", "initial_mem_free_mb", "initial_mem_used_mb",
    "post_kmer_load_mem_max_mb", "post_kmer_load_mem_total_mb",
    "post_kmer_load_mem_free_mb", "post_kmer_load_mem_used_mb",
    "input_mode",
]

CONFIG_ORDER = [
    "param_in1", "param_in2", "param_out", "param_ref", "param_stats",
    "param_k", "param_hdist", "param_mcf", "param_mm",
    "param_interleaved", "param_ordered", "param_removeifeitherbad", "param_threads",
    "jvm_flags", "bbduk_version",
    "set_interleaved", "set_ordered", "reset_interleaved", "reset_interleaved_reason",
]


def ordered_columns(present_keys, preferred_order):
    ordered = [k for k in preferred_order if k in present_keys]
    extras = [k for k in present_keys if k not in preferred_order]
    return ordered + extras


def companion_path(output_path):
    """foo.csv -> foo.run_config.csv (used for the run-parameters side file)."""
    base, ext = os.path.splitext(output_path)
    if not ext:
        base, ext = output_path, ".csv"
    return f"{base}.run_config{ext}"


def slugify(label):
    return re.sub(r"[^a-z0-9]+", "_", label.strip().lower()).strip("_")


def derive_sample_id(fields, filename):
    """Best-effort sample name: from in1= path, else the log filename stem."""
    in1 = fields.get("param_in1")
    if in1:
        base = os.path.basename(in1)
        base = re.sub(r"(_R?1)?\.(fastq|fq)(\.gz)?$", "", base, flags=re.IGNORECASE)
        if base:
            return base
    return os.path.splitext(os.path.basename(filename))[0]


def parse_command_list(text, fields):
    m = RE_COMMAND_LIST.search(text)
    if not m:
        return
    jvm_flags = []
    for token in m.group(1).split(", "):
        token = token.strip()
        if not token:
            continue
        if "=" in token:
            key, _, value = token.partition("=")
            fields[f"param_{key}"] = value
        else:
            jvm_flags.append(token)
    if jvm_flags:
        fields["jvm_flags"] = " ".join(jvm_flags)


def parse_bbduk_log(text, filename):
    fields = {"log_file": filename}
    parse_command_list(text, fields)

    m = RE_VERSION.search(text)
    if m:
        fields["bbduk_version"] = m.group(1)

    for key, val in RE_SET.findall(text):
        fields[f"set_{key.lower()}"] = val
    for key, val, reason in RE_RESET.findall(text):
        fields[f"reset_{key.lower()}"] = val
        if reason:
            fields[f"reset_{key.lower()}_reason"] = reason

    m = RE_STARTUP_TIME.search(text)
    if m:
        fields["startup_time_seconds"] = m.group(1)

    mem_blocks = RE_MEMORY.findall(text)
    if len(mem_blocks) >= 1:
        fields["initial_mem_max_mb"], fields["initial_mem_total_mb"], \
            fields["initial_mem_free_mb"], fields["initial_mem_used_mb"] = mem_blocks[0]
    if len(mem_blocks) >= 2:
        fields["post_kmer_load_mem_max_mb"], fields["post_kmer_load_mem_total_mb"], \
            fields["post_kmer_load_mem_free_mb"], fields["post_kmer_load_mem_used_mb"] = mem_blocks[1]

    m = RE_KMERS_ADDED.search(text)
    if m:
        fields["kmers_added"] = m.group(1)
        fields["kmer_load_time_seconds"] = m.group(2)

    m = RE_INPUT_MODE.search(text)
    if m:
        fields["input_mode"] = m.group(1)

    m = RE_OUTPUT_STREAMS.search(text)
    if m:
        fields["output_stream_start_seconds"] = m.group(1)

    m = RE_PROCESSING_TIME.search(text)
    if m:
        fields["processing_time_seconds"] = m.group(1)

    m = RE_TOTAL_TIME.search(text)
    if m:
        fields["total_time_seconds"] = m.group(1)

    for label, raw_count, rate, unit in RE_RATE_LINE.findall(text):
        key = slugify(label)
        fields[f"{key}_count"] = raw_count
        fields[f"{key}_rate"] = f"{rate} {unit}"

    seen_labels = set()
    for label, reads, reads_pct, bases, bases_pct in RE_STAT_PCT.findall(text):
        key = slugify(label)
        seen_labels.add(key)
        fields[f"{key}_reads"] = reads
        fields[f"{key}_reads_pct"] = reads_pct
        fields[f"{key}_bases"] = bases
        fields[f"{key}_bases_pct"] = bases_pct
    for label, reads, bases in RE_STAT_PLAIN.findall(text):
        key = slugify(label)
        if key in seen_labels:
            continue
        fields[f"{key}_reads"] = reads
        fields[f"{key}_bases"] = bases

    fields["sample_id"] = derive_sample_id(fields, filename)
    return fields


def find_candidate_files(directory, pattern, recursive):
    matches = []
    if recursive:
        for root, _dirs, files in os.walk(directory):
            for name in files:
                if fnmatch(name, pattern):
                    matches.append(os.path.join(root, name))
    else:
        for name in os.listdir(directory):
            path = os.path.join(directory, name)
            if os.path.isfile(path) and fnmatch(name, pattern):
                matches.append(path)
    return sorted(matches)


def looks_like_bbduk_log(path):
    try:
        with open(path, "r", errors="ignore") as fh:
            chunk = fh.read(20000)
        return SIGNATURE in chunk
    except OSError:
        return False


def collect_candidates(paths, pattern, recursive):
    """paths may be directories and/or individual files (e.g. a shell-expanded
    'mydir/*' glob hands us a long list of files instead of one directory)."""
    candidates = []
    seen = set()
    any_valid = False
    for p in paths:
        if os.path.isdir(p):
            any_valid = True
            for f in find_candidate_files(p, pattern, recursive=recursive):
                if f not in seen:
                    seen.add(f)
                    candidates.append(f)
        elif os.path.isfile(p):
            any_valid = True
            if fnmatch(os.path.basename(p), pattern) and p not in seen:
                seen.add(p)
                candidates.append(p)
        else:
            print(f"Warning: path not found, skipping: {p}", file=sys.stderr)
    if not any_valid:
        sys.exit("None of the given paths exist.")
    return candidates


def main():
    ap = argparse.ArgumentParser(
        description="Extract BBDuk log metadata from every matching file into one CSV.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("paths", nargs="+",
                     help="One or more directories and/or files to scan. "
                          "Shell globs like 'mydir/*' are fine -- each expanded "
                          "file is checked individually.")
    ap.add_argument("-o", "--output", default="bbduk_metadata_summary.csv",
                     help="Output CSV path (default: %(default)s)")
    ap.add_argument("-p", "--pattern", default="*",
                     help="Filename glob to consider before content-checking (default: %(default)s)")
    ap.add_argument("--no-recursive", action="store_true",
                     help="When a directory is given, only scan its top level, don't recurse")
    ap.add_argument("--single-file", action="store_true",
                     help="Write one CSV with every column instead of splitting into a "
                          "metrics matrix + a run-config reference file")
    args = ap.parse_args()

    candidates = collect_candidates(args.paths, args.pattern, recursive=not args.no_recursive)
    log_files = [p for p in candidates if looks_like_bbduk_log(p)]

    if not log_files:
        sys.exit(
            f"No BBDuk log files found among the given path(s) "
            f"(scanned {len(candidates)} file(s) matching '{args.pattern}')."
        )

    rows = [parse_bbduk_log(open(p, "r", errors="ignore").read(), p) for p in log_files]

    discovered = []
    for row in rows:
        for key in row:
            if key not in ("sample_id", "log_file") and key not in discovered:
                discovered.append(key)

    metric_keys = ordered_columns([k for k in discovered if not is_config_column(k)], METRIC_ORDER)
    config_keys = ordered_columns([k for k in discovered if is_config_column(k)], CONFIG_ORDER)

    def write_csv(path, fieldnames):
        with open(path, "w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=fieldnames, restval="", extrasaction="ignore")
            writer.writeheader()
            writer.writerows(rows)

    print(f"Parsed {len(log_files)} BBDuk log file(s).")

    if args.single_file:
        fieldnames = ["sample_id", "log_file"] + metric_keys + config_keys
        write_csv(args.output, fieldnames)
        print(f"-> {args.output}  ({len(fieldnames)} columns: sample_id x everything)")
    else:
        metrics_fieldnames = ["sample_id", "log_file"] + metric_keys
        config_fieldnames = ["sample_id", "log_file"] + config_keys
        write_csv(args.output, metrics_fieldnames)
        config_out = companion_path(args.output)
        write_csv(config_out, config_fieldnames)
        print(f"Metrics matrix   -> {args.output}  ({len(metrics_fieldnames)} columns: sample_id x metric)")
        print(f"Run config       -> {config_out}  ({len(config_fieldnames)} columns: sample_id x parameter, for reference/joins)")


if __name__ == "__main__":
    main()
