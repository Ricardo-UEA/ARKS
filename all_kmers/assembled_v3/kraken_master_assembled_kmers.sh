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

source activate upset_plot

module load krakenuniq

DB="/gpfs/data/datasets/krakendb_eupathdb54"
IN="/gpfs/home/hce24xau/scratch/gen_kmers/data/all_kmers/assembled_v3/ARKS_master_kmers_1st_screen_filtered.fasta"
OUT="/gpfs/home/hce24xau/scratch/gen_kmers/data/all_kmers/assembled_v3"
REP="/gpfs/home/hce24xau/scratch/gen_kmers/data/all_kmers/assembled_v3"

mkdir -p "$OUT" "$REP"

echo "Preloading KrakenUniq database..."
krakenuniq --db "$DB" --preload

echo "Processing $IN"

krakenuniq --db "$DB" \
    --threads 24 \
    --output "${OUT}/ARKS_1st_assembly_kraken_output" \
    --report-file "${REP}/ARKS_1st_assembly_kraken_report" \
    "$IN"

echo "Done."
