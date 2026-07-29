#!/bin/bash
#SBATCH --job-name=krakenuniq_chinese
#SBATCH --mail-type=ALL
#SBATCH --mail-user=hce24xau@uea.ac.uk
#SBATCH -o krakenuniq_chinese-%j.out
#SBATCH -e krakenuniq_chinese-%j.err
#SBATCH -p hmem-4T
#SBATCH --qos=hmem
#SBATCH --mem=500GB
#SBATCH --time=48:00:00
#SBATCH --export=ALL
#SBATCH --cpus-per-task=32

# Load environment
source activate upset_plot

# Define paths
DATABASE_PATH="/gpfs/data/datasets/krakendb_eupathdb54"
INPUT_FASTA="/gpfs/home/hce24xau/scratch/gen_kmers/data/chinese_pan_genome/output_files/chinese_pan_genome_kmers.fasta"
OUTPUT_DIR="/gpfs/home/hce24xau/scratch/gen_kmers/data/chinese_pan_genome/kraken_outputs"
REPORTS_DIR="/gpfs/home/hce24xau/scratch/gen_kmers/data/chinese_pan_genome/kraken_reports"

# Create output dirs
mkdir -p "$OUTPUT_DIR" "$REPORTS_DIR"

# Extract filename
BASENAME=$(basename "$INPUT_FASTA" .fasta)
OUTPUT_FILE="$OUTPUT_DIR/${BASENAME}_kraken_output"
REPORT_FILE="$REPORTS_DIR/${BASENAME}_kraken_report"

# Preload KrakenUniq DB
echo "Preloading KrakenUniq database..."
krakenuniq --db "$DATABASE_PATH" --preload

# Run KrakenUniq
echo "Running KrakenUniq on $INPUT_FASTA..."
krakenuniq --db "$DATABASE_PATH" \
           --threads 32 \
           --output "$OUTPUT_FILE" \
           --report-file "$REPORT_FILE" \
           "$INPUT_FASTA"

echo "Finished. Results:"
echo " - Output: $OUTPUT_FILE"
echo " - Report: $REPORT_FILE"
