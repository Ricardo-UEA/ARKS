#!/bin/bash
#SBATCH --job-name=african_batch1_krakenuniq_depleted
#SBATCH --mail-type=ALL
#SBATCH --mail-user=hce24xau@uea.ac.uk
#SBATCH -o /gpfs/home/hce24xau/scratch/gen_kmers/data/african_sydney/job_outputs/krakenuniq_mcf_05_depleted_%j.out
#SBATCH -e /gpfs/home/hce24xau/scratch/gen_kmers/data/african_sydney/job_errors/krakenuniq_mcf_05_depleted_%j.err
#SBATCH -p hmem
#SBATCH --qos=hmem
#SBATCH --mem=700GB
#SBATCH --time=168:00:00
#SBATCH --cpus-per-task=12

source activate upset_plot

BASE="/gpfs/home/hce24xau/scratch/gen_kmers/data/african_sydney"
DB="/gpfs/data/datasets/krakendb_eupathdb54"
INDIR="$BASE/output_files/hdist_0_mcf_05_t2t/clean"
OUTDIR="$BASE/kraken_outputs/t2t_mcf_05_depletion_outputs"
REPDIR="$BASE/kraken_reports/t2t_mcf_05_depletion_reports"

mkdir -p "$OUTDIR" "$REPDIR" "$BASE/job_outputs" "$BASE/job_errors"

echo "Running on host: $(hostname)"
echo "Started at: $(date)"
echo "Input directory: $INDIR"
echo "Output directory: $OUTDIR"
echo "Report directory: $REPDIR"
echo "Database: $DB"
echo

echo "Checking input files..."
ls "$INDIR"/*_clean.fastq.gz | head
echo

echo "Preloading database..."
krakenuniq --db "$DB" --preload

echo "Starting classifications..."

for FILE in "$INDIR"/*_clean.fastq.gz; do
    [ -e "$FILE" ] || continue

    SAMPLE=$(basename "$FILE" _clean.fastq.gz)
    OUT="$OUTDIR/${SAMPLE}.krakenuniq.out"
    REP="$REPDIR/${SAMPLE}.krakenuniq.report"

    if [ -s "$OUT" ] && [ -s "$REP" ]; then
        echo "Skipping $SAMPLE (already done)"
        continue
    fi

    echo "Processing $SAMPLE at $(date)"
    echo "Input file: $FILE"

    krakenuniq \
        --db "$DB" \
        --threads "$SLURM_CPUS_PER_TASK" \
        --fastq-input \
        --gzip-compressed \
        --output "$OUT" \
        --report-file "$REP" \
        "$FILE"

    echo "Finished $SAMPLE at $(date)"
    echo
done

echo "All done at $(date)"
