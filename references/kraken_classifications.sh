#!/bin/bash
#SBATCH --job-name=krakenuniq_refs
#SBATCH --mail-type=ALL
#SBATCH --mail-user=hce24xau@uea.ac.uk
#SBATCH -o krakenuniq_refs-%A_%a.out
#SBATCH -e krakenuniq_refs-%A_%a.err
#SBATCH -p hmem
#SBATCH --qos=hmem
#SBATCH --mem=800GB
#SBATCH --time=48:00:00
#SBATCH --export=ALL
#SBATCH --cpus-per-task=20
#SBATCH --array=0-2

module load krakenuniq

DATABASE_PATH="/gpfs/data/datasets/krakendb_eupathdb54"

FASTAS=(
"/gpfs/home/hce24xau/scratch/gen_kmers/data/references/input_files/GRCh37.fasta"
"/gpfs/home/hce24xau/scratch/gen_kmers/data/references/input_files/GRCh38.fasta"
"/gpfs/home/hce24xau/scratch/gen_kmers/data/references/input_files/T2T.fasta"
)

INPUT_FASTA=${FASTAS[$SLURM_ARRAY_TASK_ID]}

OUTPUT_DIR="/gpfs/home/hce24xau/scratch/gen_kmers/data/references/kraken_outputs"
REPORTS_DIR="/gpfs/home/hce24xau/scratch/gen_kmers/data/references/kraken_reports"

mkdir -p "$OUTPUT_DIR" "$REPORTS_DIR"

BASENAME=$(basename "$INPUT_FASTA" .fasta)

OUTPUT_FILE="$OUTPUT_DIR/${BASENAME}_kraken_output"
REPORT_FILE="$REPORTS_DIR/${BASENAME}_kraken_report"

echo "Running KrakenUniq on $BASENAME"

krakenuniq \
    --db "$DATABASE_PATH" \
    --preload \
    --threads 20 \
    --output "$OUTPUT_FILE" \
    --report-file "$REPORT_FILE" \
    "$INPUT_FASTA"

echo "Finished $BASENAME"
