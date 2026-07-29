#!/bin/bash
#SBATCH --job-name=jellyfish_unitigs
#SBATCH -o jellyfish_unitigs_%j.out
#SBATCH -e jellyfish_unitigs_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=hce24xau@uea.ac.uk
#SBATCH -p hmem
#SBATCH --qos=hmem
#SBATCH --mem=2000G
#SBATCH --cpus-per-task=24
#SBATCH --time=168:00:00
#SBATCH --export=ALL

source ~/.bashrc
conda activate upset_plot

K=31
THREADS=24
HASH=100G

INPUT="/gpfs/home/hce24xau/scratch/gen_kmers/data/all_kmers/assembled_v3/ARKS_1st_assembly.fasta"

JF="ARKS_master_1st_assembly.jf"
STATS="ARKS_master_1st_assembly.stats.txt"
DUMP="ARKS_master_1st_assembly_dump.fa"

echo "Counting distinct ${K}-mers..."

jellyfish count \
    -C
    -m $K \
    -s $HASH \
    -t $THREADS \
    "$INPUT" \
    -o "$JF"

cd /gpfs/home/hce24xau/scratch/gen_kmers/data/all_kmers/assembled_v3

mv ARKS_master_final.jf_0 > ARKS_master_final.jf

echo "Generating statistics..."
jellyfish stats "$JF" > "$STATS"

echo "Dumping kmers..."
jellyfish dump \
    -c \
    ARKS_master_1st_assembly.jf > "$DUMP"

echo "Done."
