#!/bin/bash
#SBATCH --job-name=jellyfish_refs
#SBATCH --mail-type=ALL
#SBATCH --mail-user=hce24xau@uea.ac.uk
#SBATCH -o jellyfish_refs-%j.out
#SBATCH -e jellyfish_refs-%j.err
#SBATCH -p compute
#SBATCH --mem=400GB
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=10
#SBATCH --export=ALL

BASE="/gpfs/home/hce24xau/scratch/gen_kmers/data/references"
INPUT_DIR="$BASE/input_files"
KMER_DIR="$BASE/kmers_files"

mkdir -p "$KMER_DIR"

source ~/.bashrc
conda activate upset_plot


# GRCh37
#jellyfish count -m 31 -s 30G -t 10 -C -o "$KMER_DIR/GRCh37.jf" "$INPUT_DIR/GRCh37.fasta"
jellyfish dump -c "$KMER_DIR/GRCh37.jf" | cut -f1 -d' ' > "$KMER_DIR/GRCh37.kmers.txt"

# GRCh38
#jellyfish count -m 31 -s 30G -t 10 -C -o "$KMER_DIR/GRCh38.jf" "$INPUT_DIR/GRCh38.fasta"
jellyfish dump -c "$KMER_DIR/GRCh38.jf" | cut -f1 -d' ' > "$KMER_DIR/GRCh38.kmers.txt"

# T2T
#jellyfish count -m 31 -s 30G -t 10 -C -o "$KMER_DIR/T2T.jf" "$INPUT_DIR/T2T.fasta"
jellyfish dump -c "$KMER_DIR/T2T.jf" | cut -f1 -d' ' > "$KMER_DIR/T2T.kmers.txt"

echo "Done"
