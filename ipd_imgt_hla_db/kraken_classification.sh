#!/bin/bash
#SBATCH --job-name=krakenuniq_hla
#SBATCH --mail-type=ALL
#SBATCH --mail-user=hce24xau@uea.ac.uk
#SBATCH -o krakenuniq_hla-%j.out
#SBATCH -e krakenuniq_hla-%j.err
#SBATCH -p hmem-4T
#SBATCH --qos=hmem
#SBATCH --mem=500GB
#SBATCH --time=48:00:00
#SBATCH --export=ALL
#SBATCH --cpus-per-task=32

source activate upset_plot

# Define input/output
DATABASE_PATH="/gpfs/data/datasets/krakendb_eupathdb54"
INPUT_FASTA="/gpfs/home/hce24xau/scratch/gen_kmers/data/ipd_imgt_hla_db/output_files/hla_nuc_kmers.fasta"
OUTPUT_DIR="/gpfs/home/hce24xau/scratch/gen_kmers/data/ipd_imgt_hla_db/kraken_outputs"
REPORTS_DIR="/gpfs/home/hce24xau/scratch/gen_kmers/data/ipd_imgt_hla_db/kraken_reports"
mkdir -p "$OUTPUT_DIR" "$REPORTS_DIR"

# Extract base name
BASENAME=$(basename "$INPUT_FASTA" .fasta)
OUTPUT_FILE="$OUTPUT_DIR/${BASENAME}_kraken_output"
REPORT_FILE="$REPORTS_DIR/${BASENAME}_kraken_report"

# Preload DB
echo "Preloading KrakenUniq database..."
krakenuniq --db "$DATABASE_PATH" --preload

# Run KrakenUniq
echo "Classifying $INPUT_FASTA..."
krakenuniq --db "$DATABASE_PATH" \
           --threads 32 \
           --output "$OUTPUT_FILE" \
           --report-file "$REPORT_FILE" \
           "$INPUT_FASTA"

echo "Done. Output saved to:"
echo "  $OUTPUT_FILE"
echo "  $REPORT_FILE"
