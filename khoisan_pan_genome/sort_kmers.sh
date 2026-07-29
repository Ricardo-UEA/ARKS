#!/bin/bash
#SBATCH --job-name=sort_kmers
#SBATCH --mail-type=ALL
#SBATCH --mail-user=hce24xau@uea.ac.uk
#SBATCH -o sort_kmers-%j.out
#SBATCH -e sort_kmers-%j.err
#SBATCH -p hmem
#SBATCH --qos=hmem
#SBATCH --mem=100GB
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=24
#SBATCH --export=ALL

set -euo pipefail

INPUT="/gpfs/home/hce24xau/scratch/gen_kmers/data/khoisan_pan_genome/input_files/KSPG.kmers.txt"
OUTPUT="/gpfs/home/hce24xau/scratch/gen_kmers/data/khoisan_pan_genome/input_files/KSPG.sorted.txt"
TMPDIR="/gpfs/home/hce24xau/scratch/gen_kmers/data/khoisan_pan_genome/sort_tmp"

echo "=== Sort started: $(date) ==="
echo "Input size: $(du -sh "$INPUT" | cut -f1)"

# Create temp directory for sort chunks
mkdir -p "$TMPDIR"

# GNU sort with:
# --parallel: use all available cores
# -T: temp directory for chunks (needs to be on scratch, not /tmp)
# -S: buffer size per thread
# --compress-program: compress temp chunks to save scratch space
sort \
    --parallel=24 \
    -T "$TMPDIR" \
    -S 4G \
    --compress-program=pigz \
    "$INPUT" \
    -o "$OUTPUT"

echo "=== Sort finished: $(date) ==="

# Verify line counts match
INPUT_LINES=$(wc -l < "$INPUT")
OUTPUT_LINES=$(wc -l < "$OUTPUT")
echo "Input lines:  $INPUT_LINES"
echo "Output lines: $OUTPUT_LINES"

if [ "$INPUT_LINES" -eq "$OUTPUT_LINES" ]; then
    echo "SUCCESS: Line counts match"
    # Clean up temp directory
    rm -rf "$TMPDIR"
else
    echo "WARNING: Line counts differ - do not delete input"
    exit 1
fi
