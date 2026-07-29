#!/bin/bash
#SBATCH --job-name=bam2fastq_african_batch2
#SBATCH --mail-type=ALL
#SBATCH --mail-user=hce24xau@uea.ac.uk
#SBATCH -o /gpfs/home/hce24xau/scratch/gen_kmers/data/african_sydney/job_outputs/bam2fastq/bam2fastq_%A_%a.out
#SBATCH -e /gpfs/home/hce24xau/scratch/gen_kmers/data/african_sydney/job_errors/bam2fastq/bam2fastq_%A_%a.err
#SBATCH -p compute
#SBATCH --mem=32GB
#SBATCH --time=12:00:00
#SBATCH --export=ALL
#SBATCH --cpus-per-task=6
#SBATCH --array=1-297

set -euo pipefail

module load samtools

INPUT_DIR="/gpfs/home/hce24xau/scratch/gen_kmers/data/african_sydney/input_files/african_batch2"
FASTQ_DIR="${INPUT_DIR}/fastq"
SORT_DIR="${INPUT_DIR}/namesort_bam"

mkdir -p "$FASTQ_DIR" "$SORT_DIR"
mkdir -p /gpfs/home/hce24xau/scratch/gen_kmers/data/african_sydney/job_outputs/bam2fastq
mkdir -p /gpfs/home/hce24xau/scratch/gen_kmers/data/african_sydney/job_errors/bam2fastq

cd "$INPUT_DIR"

mapfile -t BAMS < <(ls *.final.bam | sort)

INDEX=$((SLURM_ARRAY_TASK_ID - 1))

if [ "$INDEX" -ge "${#BAMS[@]}" ]; then
    echo "ERROR: Index $INDEX out of range for ${#BAMS[@]} BAM files"
    exit 1
fi

bam="${BAMS[$INDEX]}"
base="${bam%.final.bam}"

echo "SLURM_ARRAY_TASK_ID: ${SLURM_ARRAY_TASK_ID}"
echo "Processing BAM: ${bam}"
echo "Base name: ${base}"

samtools sort -n -@ "${SLURM_CPUS_PER_TASK}" -o "${SORT_DIR}/${base}.namesort.bam" "$bam"

samtools fastq \
    -@ "${SLURM_CPUS_PER_TASK}" \
    -1 "${FASTQ_DIR}/${base}_R1.fastq.gz" \
    -2 "${FASTQ_DIR}/${base}_R2.fastq.gz" \
    -0 /dev/null \
    -s /dev/null \
    -n "${SORT_DIR}/${base}.namesort.bam"

rm "${SORT_DIR}/${base}.namesort.bam"

echo "Finished ${base}"
