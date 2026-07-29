#!/bin/bash
#SBATCH --job-name=african_batch1_krakenuniq
#SBATCH --mail-type=ALL
#SBATCH --mail-user=hce24xau@uea.ac.uk
#SBATCH -o /gpfs/home/hce24xau/scratch/gen_kmers/data/african_sydney/job_outputs/krakenuniq_batch1_%j.out
#SBATCH -e /gpfs/home/hce24xau/scratch/gen_kmers/data/african_sydney/job_errors/krakenuniq/krakenuniq_batch1_%j.err
#SBATCH -p hmem
#SBATCH --qos=hmem
#SBATCH --mem=2000GB
#SBATCH --time=168:00:00
#SBATCH --cpus-per-task=12

source activate upset_plot

BASE="/gpfs/home/hce24xau/scratch/gen_kmers/data/african_sydney"
DB="/gpfs/data/datasets/krakendb_eupathdb54"
INDIR="$BASE/input_files/african_batch1"
OUTDIR="$BASE/kraken_outputs/african_batch1"
REPDIR="$BASE/kraken_reports/african_batch1"

mkdir -p "$OUTDIR" "$REPDIR" "$BASE/job_outputs" "$BASE/job_errors"

echo "Running on host: $(hostname)"
echo "Preloading database..."
krakenuniq --db "$DB" --preload

echo "Starting classifications..."

for R1 in "$INDIR"/*_R1.fastq.gz; do
    [ -e "$R1" ] || continue

    SAMPLE=$(basename "$R1" _R1.fastq.gz)
    R2="$INDIR/${SAMPLE}_R2.fastq.gz"
    OUT="$OUTDIR/${SAMPLE}.krakenuniq.out"
    REP="$REPDIR/${SAMPLE}.krakenuniq.report"

    if [ ! -f "$R2" ]; then
        echo "Missing R2 for $SAMPLE"
        continue
    fi

    if [ -s "$OUT" ] && [ -s "$REP" ]; then
        echo "Skipping $SAMPLE (already done)"
        continue
    fi

    echo "Processing $SAMPLE at $(date)"

    krakenuniq \
        --db "$DB" \
        --threads "$SLURM_CPUS_PER_TASK" \
        --paired \
        --fastq-input \
        --gzip-compressed \
        --output "$OUT" \
        --report-file "$REP" \
        "$R1" "$R2"

    echo "Finished $SAMPLE at $(date)"
    echo
done

echo "All done at $(date)"
