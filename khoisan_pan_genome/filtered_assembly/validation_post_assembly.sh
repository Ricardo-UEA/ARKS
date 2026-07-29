#!/bin/bash
#SBATCH --job-name=combine_kmers
#SBATCH --mail-type=ALL
#SBATCH --mail-user=hce24xau@uea.ac.uk
#SBATCH -o combine_kmers-%j.out
#SBATCH -e combine_kmers-%j.err
#SBATCH -p hmem-4T
#SBATCH --qos=hmem
#SBATCH --mem=900GB
#SBATCH --time=168:00:00
#SBATCH --export=ALL
#SBATCH --cpus-per-task=20

# ----------- CONFIGURATION -----------
module load jellyfish

KMER_SIZE=31
HASH_SIZE=100G
THREADS=20

OUTDIR="/gpfs/home/hce24xau/scratch/gen_kmers/data/african_pan/filtered_assembly"
OUT_PREFIX="african_pan_validation"
JF_FILE="${OUTDIR}/${OUT_PREFIX}.jf"
DUMP_FILE="${OUTDIR}/${OUT_PREFIX}_counts.txt"
STATS_FILE="${OUTDIR}/${OUT_PREFIX}_stats.txt"
HISTO_FILE="${OUTDIR}/${OUT_PREFIX}_histo.txt"

# Input (assembled FASTA)
FA1="/gpfs/home/hce24xau/scratch/gen_kmers/data/african_pan/filtered_assembly/filtered_african_pan_assembled.fasta"

mkdir -p "$OUTDIR"

# ----------- RUN PIPELINE -----------
echo "Jellyfish k-mer counting..."
jellyfish count -m "$KMER_SIZE" -s "$HASH_SIZE" -t "$THREADS" -C \
    -o "$JF_FILE" "$FA1"

echo "Dumping k-mer counts..."
jellyfish dump -c "$JF_FILE" | awk '{print $1}' > "$DUMP_FILE"

echo "Generating k-mer stats..."
jellyfish stats "$JF_FILE" > "$STATS_FILE"

echo "Generating k-mer histogram..."
jellyfish histo -t "$THREADS" "$JF_FILE" > "$HISTO_FILE"

# ----------- CLEANUP -----------
echo "Removing intermediate files..."
rm -f "$JF_FILE"

echo "Done. Final output files:"
echo "- $DUMP_FILE"
echo "- $STATS_FILE"
echo "- $HISTO_FILE"
