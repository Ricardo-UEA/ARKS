#!/bin/bash
#SBATCH --job-name=compare_kmers
#SBATCH --mail-type=ALL
#SBATCH --mail-user=hce24xau@uea.ac.uk
#SBATCH -o compare_kmers-%j.out
#SBATCH -e compare_kmers-%j.err
#SBATCH -p hmem
#SBATCH --qos=hmem
#SBATCH --mem=800GB
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=24
#SBATCH --export=ALL

set -euo pipefail

ARKS="/gpfs/home/hce24xau/scratch/gen_kmers/data/all_kmers/assembled_v3/ARKS_master_kmers_filtered_2nd_screen.txt"

AFRICAN="/gpfs/home/hce24xau/scratch/gen_kmers/data/african_sydney/input_files/african_batch1/all_samples.txt"
AFRICAN_SORTED="/gpfs/home/hce24xau/scratch/gen_kmers/data/african_sydney/input_files/african_batch1/all_samples_sorted.txt"

OUTDIR="/gpfs/home/hce24xau/scratch/gen_kmers/data/african_sydney/input_files/african_batch1/kmer_comparison"
TMPDIR="${OUTDIR}/sort_tmp"

mkdir -p "$OUTDIR"
mkdir -p "$TMPDIR"

echo "======================================"
echo "Started: $(date)"
echo "======================================"

echo "Sorting African k-mers..."
sort \
    --parallel=${SLURM_CPUS_PER_TASK} \
    -T "$TMPDIR" \
    -S 4G \
    --compress-program=pigz \
    "$AFRICAN" \
    -o "$AFRICAN_SORTED"

echo "Checking sort..."
INPUT_LINES=$(wc -l < "$AFRICAN")
OUTPUT_LINES=$(wc -l < "$AFRICAN_SORTED")

echo "Input lines : $INPUT_LINES"
echo "Sorted lines: $OUTPUT_LINES"

if [ "$INPUT_LINES" -ne "$OUTPUT_LINES" ]; then
    echo "ERROR: Sorting failed."
    exit 1
fi

rm -rf "$TMPDIR"

echo "Comparing k-mer sets..."

comm -23 "$ARKS" "$AFRICAN_SORTED" > "${OUTDIR}/ARKS_unique.txt"
comm -13 "$ARKS" "$AFRICAN_SORTED" > "${OUTDIR}/African_unique.txt"
comm -12 "$ARKS" "$AFRICAN_SORTED" > "${OUTDIR}/Shared_kmers.txt"

ARKS_TOTAL=$(wc -l < "$ARKS")
AFRICAN_TOTAL=$(wc -l < "$AFRICAN_SORTED")
ARKS_UNIQUE=$(wc -l < "${OUTDIR}/ARKS_unique.txt")
AFRICAN_UNIQUE=$(wc -l < "${OUTDIR}/African_unique.txt")
SHARED=$(wc -l < "${OUTDIR}/Shared_kmers.txt")

cat <<EOF

======================================
Comparison complete
======================================
ARKS total           : $ARKS_TOTAL
African total        : $AFRICAN_TOTAL

Shared               : $SHARED
Unique to ARKS       : $ARKS_UNIQUE
Unique to African    : $AFRICAN_UNIQUE

Percent of ARKS shared      : $(awk -v s="$SHARED" -v t="$ARKS_TOTAL" 'BEGIN{printf "%.2f",100*s/t}')%
Percent of African shared   : $(awk -v s="$SHARED" -v t="$AFRICAN_TOTAL" 'BEGIN{printf "%.2f",100*s/t}')%

Output files:
${OUTDIR}/Shared_kmers.txt
${OUTDIR}/ARKS_unique.txt
${OUTDIR}/African_unique.txt
======================================

EOF

echo "Finished: $(date)"
