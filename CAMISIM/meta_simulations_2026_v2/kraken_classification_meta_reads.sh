#!/bin/bash
#SBATCH --job-name=kraken_all
#SBATCH --mail-type=ALL
#SBATCH --mail-user=hce24xau@uea.ac.uk
#SBATCH -o job_errors/kraken_all_%j.out
#SBATCH -e job_outputs/kraken_all_%j.err
#SBATCH -p hmem
#SBATCH --qos=hmem
#SBATCH --mem=700GB
#SBATCH --time=168:00:00
#SBATCH --cpus-per-task=24

# ============================================================
# KrakenUniq classification of all merged metagenome reads
# Preloads database once, then processes all samples
# ============================================================

# Directories
BASE_DIR=/gpfs/home/hce24xau/scratch/gen_kmers/data/simulations/CAMISIM/meta_simulations_2026_v2
READS_DIR=${BASE_DIR}/meta_reads
REPORTS_DIR=${BASE_DIR}/kraken_reports
OUTPUTS_DIR=${BASE_DIR}/kraken_outputs
DATABASE_PATH=/gpfs/data/datasets/krakendb_eupathdb54

# Create output directories
mkdir -p "$REPORTS_DIR"
mkdir -p "$OUTPUTS_DIR"

echo "========================================"
echo "KrakenUniq Classification - All Samples"
echo "========================================"
echo ""
echo "Reads directory: $READS_DIR"
echo "Reports directory: $REPORTS_DIR"
echo "Outputs directory: $OUTPUTS_DIR"
echo "Database: $DATABASE_PATH"
echo ""

# Activate conda environment
source activate upset_plot

# Preload database into memory
echo "Preloading KrakenUniq database..."
START_PRELOAD=$(date +%s)
krakenuniq --db "$DATABASE_PATH" --preload
END_PRELOAD=$(date +%s)
echo "Database preloaded in $((END_PRELOAD - START_PRELOAD)) seconds."
echo ""

# Track progress
TOTAL_SAMPLES=$(ls -1 ${READS_DIR}/meta_sample_*_R1.fq.gz 2>/dev/null | wc -l)
CURRENT=0
FAILED=0

echo "Found $TOTAL_SAMPLES samples to process"
echo ""

# Loop through all samples
for R1 in ${READS_DIR}/meta_sample_*_R1.fq.gz; do
    [[ -e "$R1" ]] || continue
    
    # Extract sample name
    SAMPLE=$(basename "$R1" _R1.fq.gz)
    R2="${READS_DIR}/${SAMPLE}_R2.fq.gz"
    
    CURRENT=$((CURRENT + 1))
    
    echo "----------------------------------------"
    echo "[$CURRENT/$TOTAL_SAMPLES] Processing $SAMPLE"
    echo "----------------------------------------"
    
    # Check R2 exists
    if [[ ! -f "$R2" ]]; then
        echo "ERROR: Missing R2 for $SAMPLE"
        FAILED=$((FAILED + 1))
        continue
    fi
    
    # Output files
    OUTFILE="${OUTPUTS_DIR}/${SAMPLE}.kraken"
    REPORTFILE="${REPORTS_DIR}/${SAMPLE}.report"
    
    # Skip if already processed
    if [[ -f "$REPORTFILE" ]] && [[ -s "$REPORTFILE" ]]; then
        echo "Already processed, skipping..."
        continue
    fi
    
    echo "  R1: $R1"
    echo "  R2: $R2"
    
    START_TIME=$(date +%s)
    
    krakenuniq --db "$DATABASE_PATH" \
               --threads ${SLURM_CPUS_PER_TASK} \
               --paired \
               --output "$OUTFILE" \
               --report-file "$REPORTFILE" \
               "$R1" "$R2"
    
    END_TIME=$(date +%s)
    ELAPSED=$((END_TIME - START_TIME))
    
    echo "  Completed in ${ELAPSED} seconds"
    echo ""
done

echo "========================================"
echo "ALL SAMPLES COMPLETE"
echo "========================================"
echo "Total samples: $TOTAL_SAMPLES"
echo "Failed: $FAILED"
echo ""
echo "Reports: $REPORTS_DIR"
echo "Outputs: $OUTPUTS_DIR"
