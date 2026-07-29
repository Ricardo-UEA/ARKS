#!/bin/bash
#SBATCH --job-name=process_kmers         # Job name
#SBATCH --mail-type=ALL                  # Mail notifications (ALL, BEGIN, END, FAIL)
#SBATCH --mail-user=hce24xau@uea.ac.uk   # Email address for notifications
#SBATCH -o process_kmers-%j.out          # Standard output log (%j for job ID)
#SBATCH -e process_kmers-%j.err          # Standard error log (%j for job ID)
#SBATCH -p hmem                      # Which queue to use
#SBATCH --qos=hmem
#SBATCH --mem=2000GB                      # Memory required (adjust as needed)
#SBATCH --time=168:00:00                  # Time limit (adjust based on your job's estimated runtime)
#SBATCH --export=ALL                     # Export all environment variables
#SBATCH --cpus-per-task=24   


# Load the required Python environment (or activate virtual environment)
module load python/anaconda/2020.11
source activate upset_plot

cd /gpfs/home/hce24xau/scratch/gen_kmers/data/all_kmers 

python upset_plot_final.py
