#!/bin/bash
#SBATCH --job-name=kraken_arab_pan
#SBATCH --mail-type=ALL
#SBATCH --mail-user=hce24xau@uea.ac.uk
#SBATCH -o kraken_arab_pan-%j.out
#SBATCH -e kraken_arab_pan-%j.err
#SBATCH -p hmem
#SBATCH --qos=hmem
#SBATCH --mem=1000GB
#SBATCH --time=48:00:00
#SBATCH --export=ALL
#SBATCH --cpus-per-task=24

FA="/gpfs/home/hce24xau/scratch/gen_kmers/data/arab_pan_genome/output_files/arab_pan_genome.fasta"
DATABASE_PATH="/gpfs/data/datasets/krakendb_eupathdb54"
KRAKEN_OUTPUT="/gpfs/home/hce24xau/scratch/gen_kmers/data/arab_pan_genome/kraken_outputs/arab_pan_genome_kraken_output"
KRAKEN_REPORT="/gpfs/home/hce24xau/scratch/gen_kmers/data/arab_pan_genome/kraken_reports/arab_pan_genome_kraken_report"

mkdir -p "$(dirname $KRAKEN_OUTPUT)" "$(dirname $KRAKEN_REPORT)"

module load krakenuniq
krakenuniq --db "$DATABASE_PATH" --preload
krakenuniq --db "$DATABASE_PATH" \
           --threads 24 \
           --output "$KRAKEN_OUTPUT" \
           --report-file "$KRAKEN_REPORT" \
           "$FA"

echo "Done."
