#!/bin/bash
#SBATCH --job-name=bcalm_filtered_kmers         # Job name
#SBATCH -o bcalm_filtered_%j.out                # Standard output
#SBATCH -e bcalm_filtered_%j.err                # Standard error
#SBATCH --mail-type=ALL                         # Notifications on job completion/failure
#SBATCH --mail-user=hce24xau@uea.ac.uk          # Your email
#SBATCH -p hmem                              # Use high memory partition
#SBATCH --qos=hmem
#SBATCH --mem=2000GB                            # Request 4 TB RAM
#SBATCH --cpus-per-task=24                      # Use 64 CPU cores
#SBATCH --time=168:00:00                        # Maximum time (1 week)
#SBATCH --export=ALL

module load python/anaconda/2020.11

cd /gpfs/home/hce24xau/scratch/gen_kmers/data/all_kmers/assembled_v3

python convert_to_fasta.py
