#!/bin/bash
#SBATCH --job-name=sort_refs
#SBATCH --mail-type=ALL
#SBATCH --mail-user=r.ackbersingh@uea.ac.uk
#SBATCH -o sort_refs-%A_%a.out
#SBATCH -e sort_refs-%A_%a.err
#SBATCH -p hmem
#SBATCH --qos=hmem
#SBATCH --mem=200GB
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=24
#SBATCH --array=0-2
#SBATCH --export=ALL

REFDIR="/gpfs/home/hce24xau/scratch/gen_kmers/data/references/output_files"
TMPBASE="/gpfs/home/hce24xau/scratch/gen_kmers/data/references/sort_tmp"

# Map array index to file
REFS=("GRCh37" "GRCh38" "T2T")
REF="${REFS[$SLURM_ARRAY_TASK_ID]}"

INPUT="${REFDIR}/${REF}.txt"
OUTPUT="${REFDIR}/${REF}_sorted.txt"
TMPDIR="${TMPBASE}/${REF}_tmp"

mkdir -p "${TMPDIR}"

echo "=== Sort started [${REF}]: $(date) ==="
echo "Input size: $(du -sh "$INPUT" | cut -f1)"

sort \
    --parallel=24 \
    -T "${TMPDIR}" \
    -S 4G \
    --compress-program=pigz \
    "${INPUT}" \
    -o "${OUTPUT}"

echo "=== Sort finished [${REF}]: $(date) ==="

INPUT_LINES=$(wc -l < "${INPUT}")
OUTPUT_LINES=$(wc -l < "${OUTPUT}")
echo "Input lines:  ${INPUT_LINES}"
echo "Output lines: ${OUTPUT_LINES}"

if [ "${INPUT_LINES}" -eq "${OUTPUT_LINES}" ]; then
    echo "SUCCESS: Line counts match"
    rm -rf "${TMPDIR}"
else
    echo "WARNING: Line counts differ - do not delete input"
    exit 1
fi
