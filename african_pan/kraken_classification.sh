#!/bin/bash
#SBATCH --job-name=krakenuniq_african_pan
#SBATCH --mail-type=ALL
#SBATCH --mail-user=hce24xau@uea.ac.uk
#SBATCH -o krakenuniq_african_pan-%j.out
#SBATCH -e krakenuniq_african_pan-%j.err
#SBATCH -p hmem
#SBATCH --qos=hmem
#SBATCH --mem=800GB
#SBATCH --time=48:00:00
#SBATCH --export=ALL
#SBATCH --cpus-per-task=20

# Load environment
module load krakenuniq

# Define paths
DATABASE_PATH="/gpfs/data/datasets/krakendb_eupathdb54"
INPUT_FASTA="/gpfs/home/hce24xau/scratch/gen_kmers/data/african_pan/input_files/african_pan_1.fasta"
OUTPUT_DIR="/gpfs/home/hce24xau/scratch/gen_kmers/data/african_pan/kraken_outputs"
REPORTS_DIR="/gpfs/home/hce24xau/scratch/gen_kmers/data/african_pan/kraken_reports"

# Create output directories
mkdir -p "$OUTPUT_DIR" "$REPORTS_DIR"

# Extract filename
BASENAME=$(basename "$INPUT_FASTA" .fasta)

OUTPUT_FILE="$OUTPUT_DIR/${BASENAME}_kraken_output"
REPORT_FILE="$REPORTS_DIR/${BASENAME}_kraken_report"

# Preload KrakenUniq database
echo "Preloading KrakenUniq database..."
krakenuniq --db "$DATABASE_PATH" --preload

# Run KrakenUniq
echo "Running KrakenUniq on $INPUT_FASTA..."

krakenuniq \
    --db "$DATABASE_PATH" \
    --threads 20 \
    --output "$OUTPUT_FILE" \
    --report-file "$REPORT_FILE" \
    "$INPUT_FASTA"

echo "Finished."
echo "Output : $OUTPUT_FILE"
echo "Report : $REPORT_FILE"
