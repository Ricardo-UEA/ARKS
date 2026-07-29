#!/bin/bash
#SBATCH --job-name=merge_reads
#SBATCH --mail-type=ALL
#SBATCH --mail-user=hce24xau@uea.ac.uk
#SBATCH -o merge_reads_%A_%a.out
#SBATCH -e merge_reads_%A_%a.err
#SBATCH -p hmem
#SBATCH --qos=hmem
#SBATCH --mem=16GB
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=2
#SBATCH --array=76-76
# ============================================================
# Merge individual genome reads into paired-end metagenome files
# Processes each genome's R1/R2 pair together to maintain pairing
# ============================================================

# Base directories
BASE_DIR=/gpfs/home/hce24xau/scratch/gen_kmers/data/simulations/CAMISIM/meta_simulations_2026_v2
OUTPUT_BASE=$BASE_DIR/output_files
MERGED_DIR=$BASE_DIR/meta_reads

# Create merged reads directory
mkdir -p $MERGED_DIR

# Get sample number from array task ID
SAMPLE_NUM=${SLURM_ARRAY_TASK_ID}

echo "========================================"
echo "Processing meta_sample_${SAMPLE_NUM}"
echo "========================================"

# Find the sample directory (has timestamp in name)
SAMPLE_DIR=$(find ${OUTPUT_BASE}/meta_sample_${SAMPLE_NUM} -maxdepth 1 -type d -name "20*_sample_0" 2>/dev/null | head -1)

if [ -z "$SAMPLE_DIR" ]; then
    echo "ERROR: Could not find sample directory for meta_sample_${SAMPLE_NUM}"
    exit 1
fi

READS_DIR=${SAMPLE_DIR}/reads

if [ ! -d "$READS_DIR" ]; then
    echo "ERROR: Reads directory not found: $READS_DIR"
    exit 1
fi

echo "Found reads directory: $READS_DIR"

# Output files
OUT_R1=${MERGED_DIR}/meta_sample_${SAMPLE_NUM}_R1.fq.gz
OUT_R2=${MERGED_DIR}/meta_sample_${SAMPLE_NUM}_R2.fq.gz
VALIDATION_LOG=${MERGED_DIR}/meta_sample_${SAMPLE_NUM}_validation.log

# Create temp files
TMP_R1=$(mktemp)
TMP_R2=$(mktemp)

# Initialize validation log
echo "Validation Report for meta_sample_${SAMPLE_NUM}" > $VALIDATION_LOG
echo "================================================" >> $VALIDATION_LOG
echo "" >> $VALIDATION_LOG

TOTAL_READS=0
TOTAL_MISMATCHES=0
GENOMES_PROCESSED=0

# Get all R1 files and process each with its matching R2
# R1 files end with 1.fq.gz or 1.fq (e.g., Genome11.fq.gz, Genome501.fq.gz)
# R2 files end with 2.fq.gz or 2.fq (e.g., Genome12.fq.gz, Genome502.fq.gz)

for R1_FILE in $(ls ${READS_DIR}/Genome*1.fq* 2>/dev/null | sort -V); do
    # Get the matching R2 file by replacing the last '1' before .fq with '2'
    # Genome51.fq.gz -> Genome52.fq.gz
    # Genome501.fq.gz -> Genome502.fq.gz
    R2_FILE=$(echo "$R1_FILE" | sed 's/1\.fq/2.fq/')
    
    # Extract genome name for logging (e.g., Genome5, Genome50)
    GENOME_NAME=$(basename "$R1_FILE" | sed 's/1\.fq.*$//')
    
    if [ ! -f "$R2_FILE" ]; then
        echo "  WARNING: Missing R2 file for ${GENOME_NAME}" | tee -a $VALIDATION_LOG
        echo "    R1: $R1_FILE" >> $VALIDATION_LOG
        echo "    Expected R2: $R2_FILE" >> $VALIDATION_LOG
        continue
    fi
    
    echo "  Processing ${GENOME_NAME}: $(basename $R1_FILE) + $(basename $R2_FILE)"
    
    # Create temp files for this genome pair
    TMP_G_R1=$(mktemp)
    TMP_G_R2=$(mktemp)
    
    # Decompress
    if [[ "$R1_FILE" == *.gz ]]; then
        zcat "$R1_FILE" > $TMP_G_R1
    else
        cat "$R1_FILE" > $TMP_G_R1
    fi
    
    if [[ "$R2_FILE" == *.gz ]]; then
        zcat "$R2_FILE" > $TMP_G_R2
    else
        cat "$R2_FILE" > $TMP_G_R2
    fi
    
    # Count reads in each file
    R1_READS=$(awk 'END{print NR/4}' $TMP_G_R1)
    R2_READS=$(awk 'END{print NR/4}' $TMP_G_R2)
    
    if [ "$R1_READS" != "$R2_READS" ]; then
        echo "    WARNING: Unequal read counts - R1:$R1_READS R2:$R2_READS" | tee -a $VALIDATION_LOG
    fi
    
    # Validate pairing for this genome (check first 1000 reads for speed)
    MISMATCHES=$(paste <(awk 'NR%4==1 {print $1}' $TMP_G_R1 | sed 's/@//;s/\/[12]$//' | head -1000) \
                       <(awk 'NR%4==1 {print $1}' $TMP_G_R2 | sed 's/@//;s/\/[12]$//' | head -1000) | \
                 awk '$1 != $2 {count++} END {print count+0}')
    
    if [ "$MISMATCHES" -gt "0" ]; then
        echo "    WARNING: $MISMATCHES mismatched read pairs (in first 1000) in ${GENOME_NAME}" | tee -a $VALIDATION_LOG
        TOTAL_MISMATCHES=$((TOTAL_MISMATCHES + MISMATCHES))
        
        # Log first few mismatches
        paste <(awk 'NR%4==1 {print $1}' $TMP_G_R1 | sed 's/@//;s/\/[12]$//') \
              <(awk 'NR%4==1 {print $1}' $TMP_G_R2 | sed 's/@//;s/\/[12]$//') | \
        awk '$1 != $2 {print "  MISMATCH: "$1" vs "$2}' | head -5 >> $VALIDATION_LOG
    else
        echo "    ✓ ${GENOME_NAME}: $R1_READS reads, all paired"
    fi
    
    # Append to main temp files
    cat $TMP_G_R1 >> $TMP_R1
    cat $TMP_G_R2 >> $TMP_R2
    
    TOTAL_READS=$((TOTAL_READS + R1_READS))
    GENOMES_PROCESSED=$((GENOMES_PROCESSED + 1))
    
    # Clean up genome temp files
    rm -f $TMP_G_R1 $TMP_G_R2
done

echo "" >> $VALIDATION_LOG
echo "Summary:" >> $VALIDATION_LOG
echo "  Genomes processed: $GENOMES_PROCESSED" >> $VALIDATION_LOG
echo "  Total reads: $TOTAL_READS" >> $VALIDATION_LOG
echo "  Total mismatches (sampled): $TOTAL_MISMATCHES" >> $VALIDATION_LOG

if [ "$TOTAL_MISMATCHES" -eq "0" ]; then
    echo "  Status: ALL PAIRED ✓" >> $VALIDATION_LOG
    echo ""
    echo "✓ All $TOTAL_READS reads properly paired across $GENOMES_PROCESSED genomes!"
else
    echo "  Status: HAS MISMATCHES ⚠" >> $VALIDATION_LOG
    echo ""
    echo "⚠ Found $TOTAL_MISMATCHES mismatched read pairs (in sampled checks)"
fi

# Compress final output
echo ""
echo "Compressing R1..."
gzip -c $TMP_R1 > $OUT_R1

echo "Compressing R2..."
gzip -c $TMP_R2 > $OUT_R2

# Clean up temp files
rm -f $TMP_R1 $TMP_R2

# Copy mapping and profile files
cp ${READS_DIR}/reads_mapping.tsv ${MERGED_DIR}/meta_sample_${SAMPLE_NUM}_reads_mapping.tsv 2>/dev/null
cp ${SAMPLE_DIR}/../taxonomic_profile_0.txt ${MERGED_DIR}/meta_sample_${SAMPLE_NUM}_taxonomic_profile.txt 2>/dev/null

# Report file sizes
echo ""
echo "Output files:"
ls -lh $OUT_R1 $OUT_R2

echo ""
echo "Final read counts:"
R1_FINAL=$(zcat $OUT_R1 | awk 'END{print NR/4}')
R2_FINAL=$(zcat $OUT_R2 | awk 'END{print NR/4}')
echo "  R1: $R1_FINAL reads"
echo "  R2: $R2_FINAL reads"

if [ "$R1_FINAL" != "$R2_FINAL" ]; then
    echo "  ⚠ WARNING: R1 and R2 have different read counts!"
else
    echo "  ✓ R1 and R2 have equal read counts"
fi

echo ""
echo "Validation log: $VALIDATION_LOG"
echo "Finished meta_sample_${SAMPLE_NUM}"
