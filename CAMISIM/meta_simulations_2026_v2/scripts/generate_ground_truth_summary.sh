#!/bin/bash
#SBATCH --job-name=cleanup_camisim
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=hce24xau@uea.ac.uk
#SBATCH -o job_outputs/cleanup_%j.out
#SBATCH -e job_errors/cleanup_%j.err
#SBATCH -p hmem
#SBATCH --qos hmem
#SBATCH --mem=200GB
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=12


python /gpfs/home/hce24xau/scratch/gen_kmers/data/simulations/CAMISIM/meta_simulations_2026_v2/scripts/generate_ground_truth_summary.py
