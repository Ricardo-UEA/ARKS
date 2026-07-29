#!/bin/bash
#SBATCH --job-name=combine_kmers              # Job name
#SBATCH --mail-type=ALL
#SBATCH --mail-user=hce24xau@uea.ac.uk
#SBATCH -o combine_kmers-%j.out
#SBATCH -e combine_kmers-%j.err
#SBATCH -p hmem-4T
#SBATCH --qos=hmem
#SBATCH --mem=3000GB
#SBATCH --time=168:00:00
#SBATCH --export=ALL
#SBATCH --cpus-per-task=64

# ----------- CONFIGURATION -----------
module load jellyfish
module load python/anaconda/2019.10


KMER_SIZE=31
HASH_SIZE=100G
THREADS=64
OUTDIR="/gpfs/home/hce24xau/scratch/gen_kmers/data/all_kmers/output_files"
OUT_PREFIX="merged_kmers_july_with_cosmic"
JF_FILE="${OUTDIR}/${OUT_PREFIX}.jf"
DUMP_FILE="${OUTDIR}/${OUT_PREFIX}_counts.txt"
STATS_FILE="${OUTDIR}/${OUT_PREFIX}_stats.txt"

# Input files
FA1="/gpfs/home/hce24xau/scratch/gen_kmers/data/all_kmers/assembled/filtered_kmers_with_cosmic.fa"
FA2="/gpfs/home/hce24xau/scratch/gen_kmers/data/ipd_imgt_hla_db/input_files/hla_nuc.fasta"

# ----------- RUN PIPELINE -----------

echo "Jellyfish k-mer counting..."
jellyfish count -m $KMER_SIZE -s $HASH_SIZE -t $THREADS \
    -o $JF_FILE $FA1 $FA2

echo "Dumping k-mer counts..."
jellyfish dump -c $JF_FILE > $DUMP_FILE

echo "Generating k-mer stats..."
jellyfish stats $JF_FILE > $STATS_FILE

# ----------- CLEANUP -----------

echo "Removing intermediate files..."
rm -f $JF_FILE

echo "Done. Final output files:"
echo "- $DUMP_FILE"
echo "- $STATS_FILE"
