#!/bin/bash
#SBATCH --job-name=kmer_kraken_arab_pan
#SBATCH --mail-type=ALL
#SBATCH --mail-user=hce24xau@uea.ac.uk
#SBATCH -o kmer_kraken_arab_pan-%j.out
#SBATCH -e kmer_kraken_arab_pan-%j.err
#SBATCH -p hmem
#SBATCH --qos=hmem
#SBATCH --mem=20GB
#SBATCH --time=168:00:00
#SBATCH --export=ALL
#SBATCH --cpus-per-task=6


cd /gpfs/home/hce24xau/scratch/gen_kmers/data/khoisan_pan_genome/output_files

python convert_to_fasta.py
