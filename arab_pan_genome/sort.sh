#!/bin/bash
#SBATCH --job-name=sort_arab_pan
#SBATCH --mail-type=ALL
#SBATCH --mail-user=hce24xau@uea.ac.uk
#SBATCH -o sort_arab_pan-%j.out
#SBATCH -e sort_arab_pan-%j.err
#SBATCH -p hmem
#SBATCH --qos=hmem
#SBATCH --mem=800GB
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=24
#SBATCH --export=ALL

set -euo pipefail

INPUT="/gpfs/home/hce24xau/scratch/gen_kmers/data/arab_pan_genome/output_files/arab_pan_genome.txt"
OUTPUT="/gpfs/home/hce24xau/scratch/gen_kmers/data/arab_pan_genome/output_files/arab_pan_genome_sorted.txt"
TMPDIR="/gpfs/home/hce24xau/scratch/gen_kmers/data/arab_pan_genome/sort_tmp"

echo "=== Sort started: $(date) ==="
echo "Input size: $(du -sh "$INPUT" | cut -f1)"

mkdir -p "$TMPDIR"

sort \
    --parallel=24 \
    -T "$TMPDIR" \
    -S 4G \
    --compress-program=pigz \
    "$INPUT" \
    -o "$OUTPUT"

echo "=== Sort finished: $(date) ==="

INPUT_LINES=$(wc -l < "$INPUT")
OUTPUT_LINES=$(wc -l < "$OUTPUT")
echo "Input lines:  $INPUT_LINES"
echo "Output lines: $OUTPUT_LINES"

if [ "$INPUT_LINES" -eq "$OUTPUT_LINES" ]; then
    echo "SUCCESS: Line counts match"
    rm -rf "$TMPDIR"
else
    echo "WARNING: Line counts differ - do not delete input"
    exit 1
fi
