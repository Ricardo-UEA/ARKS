#!/bin/bash
#SBATCH --job-name=bbduk_african_reads_panhuman
#SBATCH --mail-type=ALL
#SBATCH --mail-user=hce24xau@uea.ac.uk
#SBATCH -o /gpfs/home/hce24xau/scratch/gen_kmers/data/african_sydney/jellyfish_%A_%a.out
#SBATCH -e /gpfs/home/hce24xau/scratch/gen_kmers/data/african_sydney/jellyfish_%A_%a.err
#SBATCH -p hmem
#SBATCH --qos=hmem
#SBATCH --mem=2000GB
#SBATCH --time=24:00:00
#SBATCH --export=ALL
#SBATCH --cpus-per-task=24

source activate upset_plot

cd /gpfs/home/hce24xau/scratch/gen_kmers/data/african_sydney/input_files/african_batch1

# Count 31-mers from all FASTQ files
#jellyfish count \
#    -m 31 \
#    -s 100G \
#    -t 48 \
#    -C \
#    -o all_samples.jf \
#    <(zcat *.fastq.gz)

# Dump k-mer counts
jellyfish dump \
    -c \
    all_samples.jf \
    > all_samples.dump.txt

# (Optional) Keep only the k-mer sequences (remove counts)
cut -d' ' -f1 all_samples.dump.txt > all_samples.kmers

# (Optional) Remove intermediate files
rm all_samples.dump.txt
# rm all_samples.jf    # Uncomment if you no longer need the Jellyfish database
