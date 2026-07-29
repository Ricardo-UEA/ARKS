#!/bin/bash
#SBATCH --job-name=kmer_kraken_arab_pan
#SBATCH --mail-type=ALL
#SBATCH --mail-user=hce24xau@uea.ac.uk
#SBATCH -o kmer_kraken_arab_pan-%j.out
#SBATCH -e kmer_kraken_arab_pan-%j.err
#SBATCH -p hmem
#SBATCH --qos=hmem
#SBATCH --mem=2000GB
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

KMER_FASTA="${BASE}/output_files/arab_pan_genome.fasta"

KRAKEN_OUTPUT="${BASE}/kraken_outputs/${OUT_PREFIX}_kraken_output"
KRAKEN_REPORT="${BASE}/kraken_reports/${OUT_PREFIX}_kraken_report"

DATABASE_PATH="/gpfs/data/datasets/krakendb_eupathdb54"

# ----------- START -----------
echo "=========================================="
echo "Starting k-mer KrakenUniq workflow"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Host: $(hostname)"
echo "Start time: $(date)"
echo "=========================================="

# ----------- SETUP -----------
echo
echo "[1/5] Activating conda environment..."
eval "$(conda shell.bash hook)"
source activate upset_plot

echo "[INFO] Creating output directories..."
mkdir -p "$KMER_DIR" "$DUMP_DIR" "${BASE}/output_files" \
         "${BASE}/kraken_outputs" "${BASE}/kraken_reports"

echo "[INFO] Input FASTA:"
echo "$FA"

echo "[INFO] Using:"
echo "  K-mer size : $KMER_SIZE"
echo "  Hash size  : $HASH_SIZE"
echo "  Threads    : $THREADS"

# ----------- JELLYFISH COUNT -----------
echo
echo "[2/5] Running Jellyfish count..."
echo "Output JF file:"
echo "$JF_FILE"

#jellyfish count \
#    -m $KMER_SIZE \
#    -s $HASH_SIZE \
#    -t $THREADS \
#    -C \
#    -o "$JF_FILE" \
#    "$FA"

echo "[INFO] Jellyfish count complete."

# ----------- JELLYFISH DUMP -----------
echo
echo "[3/5] Dumping k-mers from Jellyfish database..."
echo "Dump file:"
echo "$DUMP_FILE"

#jellyfish dump -c "$JF_FILE" > "$DUMP_FILE"

echo "[INFO] Removing temporary Jellyfish database..."
rm -f "$JF_FILE"

# ----------- CONVERT TO FASTA -----------
echo
echo "[4/5] Converting dumped k-mers to FASTA format..."
echo "FASTA output:"
echo "$KMER_FASTA"

#awk '{print ">kmer_" NR "_count_" $2 "\n" $1}' "$DUMP_FILE" > "$KMER_FASTA"

echo "[INFO] FASTA conversion complete."

# ----------- KRAKENUNIQ -----------
echo
echo "[5/5] Running KrakenUniq classification..."
echo "Loading KrakenUniq module..."

module load krakenuniq

echo "[INFO] Preloading KrakenUniq database..."
krakenuniq --db "$DATABASE_PATH" --preload

echo "[INFO] Starting KrakenUniq classification..."
echo "Kraken output:"
echo "$KRAKEN_OUTPUT"

echo "Kraken report:"
echo "$KRAKEN_REPORT"

krakenuniq \
    --db "$DATABASE_PATH" \
    --threads $THREADS \
    --output "$KRAKEN_OUTPUT" \
    --report-file "$KRAKEN_REPORT" \
    "$KMER_FASTA"

# ----------- COMPLETE -----------
echo
echo "=========================================="
echo "Workflow complete."
echo "Finished at: $(date)"
echo "=========================================="
