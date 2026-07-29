#!/bin/bash
#SBATCH --job-name=jellyfish_kspg
#SBATCH -o jellyfish_kspg_%j.out
#SBATCH -e jellyfish_kspg_%j.err
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

INPUT="/gpfs/home/hce24xau/scratch/gen_kmers/data/khoisan_pan_genome/output_files/KSPG.fasta"
OUTDIR="/gpfs/home/hce24xau/scratch/gen_kmers/data/khoisan_pan_genome/output_files"

JF="${OUTDIR}/KSPG.jf"
STATS="${OUTDIR}/KSPG.stats.txt"
DUMP="${OUTDIR}/KSPG.dump.txt"

echo "Counting distinct ${K}-mers..."

jellyfish count \
    -m $K \
    -s $HASH \
    -C \
    -t $THREADS \
    "$INPUT" \
    -o "$JF"

# Jellyfish may create split files
if [ -f "${JF}_0" ]; then
    mv "${JF}_0" "$JF"
fi

echo "Generating statistics..."
jellyfish stats "$JF" > "$STATS"

echo "Dumping kmers..."
jellyfish dump -c "$JF" > "$DUMP"

echo "Done."
