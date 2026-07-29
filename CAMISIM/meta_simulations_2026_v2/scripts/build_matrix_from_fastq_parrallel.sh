#!/bin/bash
#SBATCH --job-name=cleanup_camisim
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=hce24xau@uea.ac.uk
#SBATCH -o cleanup_%j.out
#SBATCH -e cleanup_%j.err
#SBATCH -p hmem
#SBATCH --qos hmem
#SBATCH --mem=200GB
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=16

python /gpfs/home/hce24xau/scratch/gen_kmers/data/simulations/CAMISIM/meta_simulations_2026_v2/scripts/build_matrix_from_fastq_parrallel.py
