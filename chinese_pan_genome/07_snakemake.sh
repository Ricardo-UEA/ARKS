#!/bin/bash
#SBATCH --job-name=process_kmers              # Job name
#SBATCH --mail-type=ALL                       # Mail notifications (ALL, BEGIN, END, FAIL)
#SBATCH --mail-user=hce24xau@uea.ac.uk        # Email address for notifications
#SBATCH -o process_kmers-%j.out               # Standard output log (%j for job ID)
#SBATCH -e process_kmers-%j.err               # Standard error log (%j for job ID)
#SBATCH -p hmem-4T 
#SBATCH --qos=hmem
#SBATCH --mem=20GB                           # Request 700GB of RAM
#SBATCH --time=168:00:00                       # Maximum time limit (72 hours)
#SBATCH --export=ALL                          # Export environment variables
#SBATCH --cpus-per-task=5                   # 24 CPU cores per task

module load jellyfish/2.3.0
module load Snakemake
module load python/anaconda/2019.10/3.7

# Define directories (make sure these paths match your actual workflow)
SNAKE_SCRIPT="/gpfs/home/hce24xau/scratch/gen_kmers/data/chinese_pan_genome/07_snakemake.smk"
PROFILE="/gpfs/home/hce24xau/.config/snakemake/slurm/profile"

# Call Snakemake with cluster submission
snakemake -s $SNAKE_SCRIPT --profile $PROFILE  --jobs 100 --latency-wait 10 
