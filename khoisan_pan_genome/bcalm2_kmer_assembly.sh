#!/bin/bash
#SBATCH --job-name=bcalm_filtered_kmers         # Job name
#SBATCH -o bcalm_filtered_%j.out                # Standard output
#SBATCH -e bcalm_filtered_%j.err                # Standard error
#SBATCH --mail-type=ALL                         # Notifications on job completion/failure
#SBATCH --mail-user=hce24xau@uea.ac.uk          # Your email
#SBATCH -p hmem                              # Use high memory partition
#SBATCH --qos=hmem
#SBATCH --mem=2000GB                            # Request 4 TB RAM
#SBATCH --cpus-per-task=24                      # Use 64 CPU cores
#SBATCH --time=168:00:00                        # Maximum time (1 week)
#SBATCH --export=ALL

echo "Starting BCALM on filtered k-mers..."

# Load conda and activate environment properly
source ~/.bashrc
conda activate bcalm_env

# Define paths
INPUT_KMERS="/gpfs/home/hce24xau/scratch/gen_kmers/data/khoisan_pan_genome/output_files/KPSG.fasta"
OUTPUT_PREFIX="/gpfs/home/hce24xau/scratch/gen_kmers/data/khoisan_pan_genome/output_files/KPSG_assembly.fasta"
TMP_DIR="/gpfs/home/hce24xau/scratch/gen_kmers/data/all_kmers/assembled_v3/tmp"

# Ensure output and temp directories exist
mkdir -p "$(dirname "$OUTPUT_PREFIX")"
mkdir -p "$TMP_DIR"

# Run BCALM
bcalm \
  -in "$INPUT_KMERS" \
  -kmer-size 31 \
  -abundance-min 1 \
  -out "$OUTPUT_PREFIX" \
  -nb-cores 24 \
  -debloom none \
  -max-memory 2000000 \
  -out-tmp "$TMP_DIR" \
  -verbose 2

echo "BCALM assembly complete."
