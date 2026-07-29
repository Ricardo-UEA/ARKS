#!/bin/bash
#SBATCH --job-name=african_pan_kraken
#SBATCH --mail-type=ALL
#SBATCH --mail-user=hce24xau@uea.ac.uk
#SBATCH -o master_pan_kraken-%j.out
#SBATCH -e master_pan_kraken-%j.err
#SBATCH -p hmem
#SBATCH --qos=hmem
#SBATCH --mem=2000GB
#SBATCH --time=168:00:00
#SBATCH --cpus-per-task=24

module load python/anaconda/2020.11

cd /gpfs/home/hce24xau/scratch/gen_kmers/data/all_kmers/assembled_v3

python extract_taxa_kmers.py \
  --report classified_kmers_kraken.tsv \
  --fasta /gpfs/home/hce24xau/scratch/gen_kmers/data/all_kmers/assembled_v3/ARKS_master_kmers_1st_screen_filtered.fasta \
  --outdir /gpfs/home/hce24xau/scratch/gen_kmers/data/all_kmers/assembled_v3/taxa_kmers_to_remove

