#!/bin/bash
#SBATCH --job-name=kmer_kraken_arab_pan
#SBATCH --mail-type=ALL
#SBATCH --mail-user=hce24xau@uea.ac.uk
#SBATCH -o kmer_kraken_arab_pan-%j.out
#SBATCH -e kmer_kraken_arab_pan-%j.err
#SBATCH -p hmem
#SBATCH --qos=hmem
#SBATCH --mem=1000GB
#SBATCH --time=168:00:00
#SBATCH --export=ALL
#SBATCH --cpus-per-task=24

# ----------- CONFIGURATION -----------
KMER_SIZE=31
HASH_SIZE=50G
THREADS=24

BASE="/gpfs/home/hce24xau/scratch/gen_kmers/data/arab_pan_genome"
FA="${BASE}/input_files/arab_pan_genome.fasta"
KMER_DIR="${BASE}/kmers_files"
DUMP_DIR="${BASE}/dump_files"
OUT_PREFIX="arab_pan_genome"
JF_FILE="${KMER_DIR}/${OUT_PREFIX}.comb.jf"
DUMP_FILE="${DUMP_DIR}/${OUT_PREFIX}_counts.txt"
KMER_FASTA="${BASE}/output_files/arab_pan_genome_kmers.fasta"
KRAKEN_OUTPUT="${BASE}/kraken_outputs/${OUT_PREFIX}_kraken_output"
KRAKEN_REPORT="${BASE}/kraken_reports/${OUT_PREFIX}_kraken_report"
DATABASE_PATH="/gpfs/data/datasets/krakendb_eupathdb54"

# ----------- SETUP -----------
eval "$(conda shell.bash hook)"
conda activate upset_plot
mkdir -p "$KMER_DIR" "$DUMP_DIR" "${BASE}/output_files" \
         "${BASE}/kraken_outputs" "${BASE}/kraken_reports"

# ----------- JELLYFISH -----------
jellyfish count -m $KMER_SIZE -s $HASH_SIZE -t $THREADS -C -o "$JF_FILE" "$FA"

jellyfish dump -c "$JF_FILE" > "$DUMP_FILE"
rm -f "$JF_FILE"

# ----------- CONVERT TO FASTA -----------
awk '{print ">kmer_" NR "_count_" $2 "\n" $1}' "$DUMP_FILE" > "$KMER_FASTA"

# ----------- KRAKENUNIQ -----------
module load krakenuniq
krakenuniq --db "$DATABASE_PATH" --preload
krakenuniq --db "$DATABASE_PATH" \
           --threads $THREADS \
           --output "$KRAKEN_OUTPUT" \
           --report-file "$KRAKEN_REPORT" \
           "$KMER_FASTA"

echo "Done."
