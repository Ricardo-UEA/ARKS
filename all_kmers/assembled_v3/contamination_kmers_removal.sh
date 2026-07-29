#!/bin/bash
#SBATCH --job-name=filter_kmers
#SBATCH --mail-type=ALL
#SBATCH --mail-user=r.ackbersingh@uea.ac.uk
#SBATCH -o filter_kmers-%j.out
#SBATCH -e filter_kmers-%j.err
#SBATCH -p hmem
#SBATCH --qos=hmem
#SBATCH --mem=320GB
#SBATCH --time=168:00:00
#SBATCH --cpus-per-task=24
#SBATCH --export=ALL

set -euo pipefail

DIR="/gpfs/home/hce24xau/scratch/gen_kmers/data/all_kmers/assembled_v3"
MASTER="${DIR}/ARKS_pre_screen_master_kmers.txt"
REMOVE="${DIR}/all_taxa_kmers_to_remove.txt"
OUT="${DIR}/ARKS_master_kmers_1st_screen_filtered.txt"

echo "=== Filter started: $(date) ==="
echo "Master k-mers:    $(wc -l < "${MASTER}")"
echo "K-mers to remove: $(wc -l < "${REMOVE}")"

comm -23 "${MASTER}" "${REMOVE}" > "${OUT}"

echo "=== Filter finished: $(date) ==="
echo "Input k-mers:   $(wc -l < "${MASTER}")"
echo "Removed k-mers: $(wc -l < "${REMOVE}")"
echo "Output k-mers:  $(wc -l < "${OUT}")"
echo "Output size:    $(du -sh "${OUT}" | cut -f1)"
