#!/bin/bash
#SBATCH --job-name=merge_kmers
#SBATCH --mail-type=ALL
#SBATCH --mail-user=r.ackbersingh@uea.ac.uk
#SBATCH -o merge_kmers-%j.out
#SBATCH -e merge_kmers-%j.err
#SBATCH -p hmem
#SBATCH --qos=hmem
#SBATCH --mem=640GB
#SBATCH --time=168:00:00
#SBATCH --cpus-per-task=24
#SBATCH --export=ALL

set -euo pipefail

SCRATCH="/gpfs/home/hce24xau/scratch/gen_kmers"
TMPDIR="${SCRATCH}/data/all_kmers/assembled_v3/merge_tmp"
OUT="${SCRATCH}/data/all_kmers/assembled_v3/ARKS_pre_screen_master_kmers.txt"

mkdir -p "${TMPDIR}"

echo "=== Merge started: $(date) ==="

sort -m -u \
    --parallel=24 \
    -S 4G \
    -T "${TMPDIR}" \
    --compress-program=pigz \
    "${SCRATCH}/data/references/output_files/GRCh37.txt" \
    "${SCRATCH}/data/references/output_files/GRCh38.txt" \
    "${SCRATCH}/data/references/output_files/T2T.txt" \
    "${SCRATCH}/data/african_pan/output_files/african_pan.txt" \
    "${SCRATCH}/data/arab_pan_genome/output_files/arab_pan_genome.txt" \
    "${SCRATCH}/data/khoisan_pan_genome/output_files/KSPG.txt" \
    "${SCRATCH}/data/chinese_pan_genome/output_files/chinese_pan_genome.txt" \
    "${SCRATCH}/data/pan_genome_2/output_files/pan_genome_2.txt" \
    "${SCRATCH}/data/ipd_imgt_hla_db/output_files/hla_nuc_kmers.txt" \
    -o "${OUT}"

echo "=== Merge finished: $(date) ==="
echo "Total unique k-mers: $(wc -l < "${OUT}")"
echo "Output size: $(du -sh "${OUT}" | cut -f1)"

rm -rf "${TMPDIR}"
