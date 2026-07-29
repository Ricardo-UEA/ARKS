#!/bin/bash
#SBATCH --job-name=simulations_depletion
#SBATCH --mail-type=ALL
#SBATCH --mail-user=hce24xau@uea.ac.uk
#SBATCH -o job_outputs/bbduk_mcf_05/depletion_ARKS_array_%A_%a.out
#SBATCH -e job_errors/bbduk_mcf_05/depletion_ARKS__array_%A_%a.err
#SBATCH -p compute
#SBATCH --mem=200GB
#SBATCH --time=168:00:00
#SBATCH --cpus-per-task=6
#SBATCH --export=ALL
#SBATCH --array=1-100

module load bbmap

BASE_DIR=/gpfs/home/hce24xau/scratch/gen_kmers/data/simulations/CAMISIM/meta_simulations_2026_v2

# Reference
REF=/gpfs/afm/hpccancergenetics/Ricardo/ARKS_assembly/ARKS_assembly.fasta

# Directories - adjust as needed
INPUT_DIR=$BASE_DIR/meta_reads
OUTPUT_DIR=$BASE_DIR/depleted_reads
STATS_DIR=$BASE_DIR/stats_files
mkdir -p "${OUTPUT_DIR}" logs
mkdir -p "${STATS_DIR}" logs

# Sample name from array index
SAMPLE="meta_sample_${SLURM_ARRAY_TASK_ID}"

R1="${INPUT_DIR}/${SAMPLE}_R1.fq.gz"
R2="${INPUT_DIR}/${SAMPLE}_R2.fq.gz"

OUT_R1="${OUTPUT_DIR}/${SAMPLE}_depleted_R1.fastq.gz"
OUT_R2="${OUTPUT_DIR}/${SAMPLE}_depleted_R2.fastq.gz"
STATS="${STATS_DIR}/${SAMPLE}_bbduk_stats.txt"

echo "[$(date)] Processing ${SAMPLE}"

bbduk.sh \
    -Xmx170g \
    in1="${R1}" \
    in2="${R2}" \
    out1="${OUT_R1}" \
    out2="${OUT_R2}" \
    ref="${REF}" \
    k=31 \
    hdist=0 \
    mm=t \
    mcf=0.5 \
    removeifeitherbad=t \
    ordered=t \
    threads="${SLURM_CPUS_PER_TASK}" \
    stats="${STATS}" \
    overwrite=true 

echo "[$(date)] Done: ${SAMPLE}"
