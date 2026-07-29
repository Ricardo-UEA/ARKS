#!/bin/bash
#SBATCH --job-name=bcalm_filtered_kmers
#SBATCH -o bcalm_filtered_%j.out
#SBATCH -e bcalm_filtered_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=hce24xau@uea.ac.uk
#SBATCH -p hmem-4T
#SBATCH --qos=hmem
#SBATCH --mem=1000GB
#SBATCH --cpus-per-task=20
#SBATCH --time=168:00:00
#SBATCH --export=ALL

echo "Starting BCALM on filtered k-mers..."

# Load conda
source ~/.bashrc
conda activate bcalm_env

# Paths
INPUT_KMERS="/gpfs/home/hce24xau/scratch/gen_kmers/data/african_pan/filtered_assembly/filtered_african_pan.fasta"
OUTPUT_PREFIX="/gpfs/home/hce24xau/scratch/gen_kmers/data/african_pan/filtered_assembly/filtered_african_pan"
TMP_DIR="/gpfs/home/hce24xau/scratch/gen_kmers/data/african_pan/filtered_assembly/tmp"

mkdir -p "$(dirname "$OUTPUT_PREFIX")" "$TMP_DIR"

# Run BCALM (max-memory in MB, ~950GB here)
bcalm \
  -in "$INPUT_KMERS" \
  -kmer-size 31 \
  -abundance-min 1 \
  -out "$OUTPUT_PREFIX" \
  -nb-cores 20 \
  -max-memory 950000 \
  -out-tmp "$TMP_DIR" \
  -verbose 2

echo "BCALM assembly complete."
