#!/usr/bin/env bash
# make_deacon_csv.sh
# Usage: bash make_deacon_csv.sh [stats_dir] [out_csv]

set -euo pipefail

STATS_DIR="${1:-/gpfs/home/hce24xau/scratch/gen_kmers/data/african_sydney/stats_files}"
OUT_CSV="${2:-${STATS_DIR%/}/deacon_summaries.csv}"

mkdir -p "$STATS_DIR"

# Write CSV header
echo "sample,version,index,input,output,k,w,abs_threshold,rel_threshold,prefix_length,deplete,rename,seqs_in,seqs_out,seqs_out_proportion,seqs_removed,seqs_removed_proportion,bp_in,bp_out,bp_out_proportion,bp_removed,bp_removed_proportion,time,seqs_per_second,bp_per_second" > "$OUT_CSV"

shopt -s nullglob
for f in "$STATS_DIR"/*.summary.json; do
python3 <<EOF >> "$OUT_CSV"
import json, os, csv

with open("$f") as fh:
    j = json.load(fh)

sample = os.path.basename(j["input"]).split(".")[0]

row = [
    sample,
    j.get("version",""),
    j.get("index",""),
    j.get("input",""),
    j.get("output",""),
    j.get("k",""),
    j.get("w",""),
    j.get("abs_threshold",""),
    j.get("rel_threshold",""),
    j.get("prefix_length",""),
    j.get("deplete",""),
    j.get("rename",""),
    j.get("seqs_in",""),
    j.get("seqs_out",""),
    j.get("seqs_out_proportion",""),
    j.get("seqs_removed",""),
    j.get("seqs_removed_proportion",""),
    j.get("bp_in",""),
    j.get("bp_out",""),
    j.get("bp_out_proportion",""),
    j.get("bp_removed",""),
    j.get("bp_removed_proportion",""),
    j.get("time",""),
    j.get("seqs_per_second",""),
    j.get("bp_per_second","")
]

writer = csv.writer(open("$OUT_CSV", "a"))
writer.writerow(row)
EOF
done

echo "✅ CSV created:"
echo "$OUT_CSV"

