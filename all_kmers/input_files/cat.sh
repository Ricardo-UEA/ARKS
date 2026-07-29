#!/bin/bash
#SBATCH --job-name=combine_kmers              # Job name
#SBATCH --mail-type=ALL                       # Mail notifications (ALL, BEGIN, END, FAIL)
#SBATCH --mail-user=hce24xau@uea.ac.uk        # Email address for notifications
#SBATCH -o combine_kmers-%j.out               # Standard output log (%j for job ID)
#SBATCH -e combine_kmers-%j.err               # Standard error log (%j for job ID)
#SBATCH -p hmem-4T
#SBATCH --qos=hmem
#SBATCH --mem=20GB                           # Request 700GB of RAM
#SBATCH --time=168:00:00                       # Maximum time limit (72 hours)
#SBATCH --export=ALL                          # Export environment variables
#SBATCH --cpus-per-task=5                   # 24 CPU cores per task


cat \
/gpfs/home/hce24xau/scratch/gen_kmers/data/ipd_imgt_hla_db/input_files/hla_nuc.fasta \
/gpfs/home/hce24xau/scratch/gen_kmers/data/african_pan/input_files/african_pan_1.fasta \
/gpfs/home/hce24xau/scratch/gen_kmers/data/chinese_pan_genome/input_files/combined_chinese_pan_genome.fasta \
/gpfs/home/hce24xau/scratch/gen_kmers/data/references/input_files/T2T.fasta \
/gpfs/home/hce24xau/scratch/gen_kmers/data/references/input_files/GRCh37.fasta \
/gpfs/home/hce24xau/scratch/gen_kmers/data/references/input_files/GRCh38.fasta \
/gpfs/home/hce24xau/scratch/gen_kmers/data/pan_genome_2/input_files/pan_genome_2.fasta \
> /gpfs/home/hce24xau/scratch/gen_kmers/data/all_kmers/input_files/all_kmers.fasta
