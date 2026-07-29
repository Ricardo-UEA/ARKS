#!/bin/bash
#SBATCH --job-name=african_depletion_kmers
#SBATCH --mail-type=ALL
#SBATCH --mail-user=hce24xau@uea.ac.uk
#SBATCH -o african_depletion_kmers_%j.out    # Job array output file
#SBATCH -e african_depletion_kmers_%j.err    # Job array error file
#SBATCH -p compute
#SBATCH --mem=100GB
#SBATCH --time=168:00:00
#SBATCH --export=ALL
#SBATCH --cpus-per-task=5  # Increase CPUs to match Kraken threads

module load bbmap
module load samtools
module load Snakemake
module load python/anaconda/2019.10

#SNAKE_SCRIPT="/gpfs/home/hce24xau/scratch/gen_kmers/data/african_sydney/deplete_comparison.smk"
PROFILE="/gpfs/home/hce24xau/.config/snakemake/slurm/profile" 
 

# Call Snakemake with cluster submission
snakemake -s $SNAKE_SCRIPT --profile $PROFILE  --jobs 100 --latency-wait 10 

