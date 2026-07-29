#!/bin/bash
#SBATCH --job-name=jellyfish_unitigs             # Job name
#SBATCH -o jellyfish_unitigs_%j.out              # Standard output
#SBATCH -e jellyfish_unitigs_%j.err              # Standard error
#SBATCH --mail-type=ALL                          # Notifications
#SBATCH --mail-user=hce24xau@uea.ac.uk           # Your email
#SBATCH -p hmem                                # High memory partition
#SBATCH --qos=hmem
#SBATCH --mem=2000GB                              # 4 TB RAM
#SBATCH --cpus-per-task=24                        # 64 cores
#SBATCH --time=168:00:00                          # 7 days
#SBATCH --export=ALL

echo " ^=^t^d Loading Jellyfish..."

source ~/.bashrc
conda activate upset_plot

# Define variables
KMER_SIZE=31
THREADS=24
HASH_SIZE=100G
INPUT_FA="/gpfs/home/hce24xau/scratch/gen_kmers/data/arab_pan_genome/input_files/arab_pan_genome.fasta"
JF_OUT="/gpfs/home/hce24xau/scratch/gen_kmers/data/arab_pan_genome/kmers_files/arab_pan_genome.jf"
DUMP_LIST="/gpfs/home/hce24xau/scratch/gen_kmers/data/arab_pan_genome/output_files/arab_pan_genome_list.txt" 
DUMP_OUT="/gpfs/home/hce24xau/scratch/gen_kmers/data/arab_pan_genome/output_files/arab_pan_genome.txt" 
STATS_OUT="validation_stats.txt"

#echo " ^=^t^a Running Jellyfish count..."
#jellyfish count -m $KMER_SIZE -s $HASH_SIZE -t $THREADS -C $INPUT_FA -o $JF_OUT

echo " ^=^s^j Running Jellyfish stats..."
jellyfish stats $JF_OUT > $STATS_OUT

echo " ^=^s  Dumping k-mers to FASTA..."
jellyfish dump -c $JF_OUT > $DUMP_LIST
 
cut -d ' ' -f 1 $DUMP_LIST > $DUMP_OUT

 
cd /gpfs/home/hce24xau/scratch/gen_kmers/data/arab_pan_genome/output_files

python convert_to_fasta.py

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

echo " ^|^e Done. Output:"
echo "  - Binary counts: $JF_OUT"
echo "  - Stats:         $STATS_OUT"
echo "  - Dumped FASTA:  $DUMP_OUT"
