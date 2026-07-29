#!/bin/bash
#SBATCH --job-name=bcalm_chinese_pan         # Job name
#SBATCH -o bcalm_chinese_pan_%j.out          # Standard output
#SBATCH -e bcalm_chinese_pan_%j.err          # Standard error
#SBATCH --mail-type=ALL                      # Notifications on job completion/failure
#SBATCH --mail-user=hce24xau@uea.ac.uk       # Your email
#SBATCH -p hmem-4T                            # Use high memory partition
#SBATCH --qos=hmem
#SBATCH --mem=4000GB                          # Request 4 TB RAM
#SBATCH --cpus-per-task=64                    # Use 64 CPU cores
#SBATCH --time=168:00:00                      # Maximum time (1 week)
#SBATCH --export=ALL

echo "Starting BCALM on Chinese pan-genome k-mers..."

# Load conda and activate environment
source ~/.bashrc
conda activate bcalm_env

# Define paths
INPUT_KMERS="/gpfs/home/hce24xau/scratch/gen_kmers/data/chinese_pan_genome/output_files/chinese_pan_genome_kmers.fasta"
OUTPUT_PREFIX="/gpfs/home/hce24xau/scratch/gen_kmers/data/chinese_pan_genome/output_files/chinese_pan_genome_kmers"
TMP_DIR="/gpfs/home/hce24xau/scratch/gen_kmers/data/chinese_pan_genome/output_files/tmp"

# Ensure output and temp directories exist
mkdir -p "$(dirname "$OUTPUT_PREFIX")"
mkdir -p "$TMP_DIR"

# Run BCALM
bcalm \
  -in "$INPUT_KMERS" \
  -kmer-size 31 \
  -abundance-min 1 \
  -out "$OUTPUT_PREFIX" \
  -nb-cores 64 \
  -max-memory 4000000 \
  -out-tmp "$TMP_DIR" \
  -verbose 2

echo "BCALM assembly complete."
