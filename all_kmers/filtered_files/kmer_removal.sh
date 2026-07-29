#!/bin/bash
#SBATCH --job-name=african_pan_kraken          # Job name
#SBATCH --mail-type=ALL                        # Mail notifications (ALL, BEGIN, END, FAIL)
#SBATCH --mail-user=hce24xau@uea.ac.uk         # Email address for notifications
#SBATCH -o kmer_removal_%j.out           # Standard output log (%j for job ID)
#SBATCH -e kmer_removal_%j.err           # Standard error log (%j for job ID)
#SBATCH -p hmem-4T
#SBATCH --qos=hmem
#SBATCH --mem=2000GB                           # Request 1000GB of RAM
#SBATCH --time=168:00:00                       # Maximum time limit (168 hours)
#SBATCH --export=ALL                           # Export environment variables
#SBATCH --cpus-per-task=32                     # 64 CPU cores per task
 
module load python/anaconda/2019.10

# Run the first script
python /gpfs/home/hce24xau/scratch/gen_kmers/data/all_kmers/filtered_files/kmer_removal.py
 
