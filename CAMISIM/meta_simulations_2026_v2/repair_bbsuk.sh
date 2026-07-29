#!/bin/bash
#SBATCH --job-name=repair_pairs
#SBATCH --mail-type=ALL
#SBATCH --mail-user=hce24xau@uea.ac.uk
#SBATCH -p compute
#SBATCH --mem=64GB
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=3
#SBATCH --array=1-100
#SBATCH -o /gpfs/home/hce24xau/scratch/gen_kmers/data/simulations/CAMISIM/meta_simulations_2026_v2/job_outputs/repair_%A_%a.out
#SBATCH -e /gpfs/home/hce24xau/scratch/gen_kmers/data/simulations/CAMISIM/meta_simulations_2026_v2/job_errors/repair_%A_%a.err

module load bbmap

READS_DIR="/gpfs/home/hce24xau/scratch/gen_kmers/data/simulations/CAMISIM/meta_simulations_2026_v2/meta_reads"
OUT_DIR="/gpfs/home/hce24xau/scratch/gen_kmers/data/simulations/CAMISIM/meta_simulations_2026_v2/meta_reads_repaired"
STATS_DIR="/gpfs/home/hce24xau/scratch/gen_kmers/data/simulations/CAMISIM/meta_simulations_2026_v2/stats_files"

mkdir -p "$OUT_DIR" "$STATS_DIR"

ID="${SLURM_ARRAY_TASK_ID}"
SAMPLE="meta_sample_${ID}"

IN1="${READS_DIR}/${SAMPLE}_R1.fq.gz"
IN2="${READS_DIR}/${SAMPLE}_R2.fq.gz"

OUT1="${OUT_DIR}/${SAMPLE}_R1.repaired.fq.gz"
OUT2="${OUT_DIR}/${SAMPLE}_R2.repaired.fq.gz"
SING="${OUT_DIR}/${SAMPLE}.singletons.fq.gz"
STAT="${STATS_DIR}/${SAMPLE}.tsv"

echo "[$(date)] Repairing ${SAMPLE}"
[[ -s "$IN1" && -s "$IN2" ]] || { echo "Missing inputs"; exit 1; }

# capture repair.sh output to parse stats
LOG="${STATS_DIR}/${SAMPLE}.repair.log"

repair.sh rp in1="$IN1" in2="$IN2" out1="$OUT1" out2="$OUT2" outs="$SING" 2>&1 | tee "$LOG"

# minimal stats extraction
pairs=$(grep -E '^Pairs:' "$LOG" | awk '{print $2}')
singletons=$(grep -E '^Singletons:' "$LOG" | awk '{print $2}')
time_s=$(grep -E '^Time:' "$LOG" | awk '{print $(NF-1)}')  # seconds.

{
  echo -e "sample\tpairs_reads\tsingletons_reads\ttime_s\tout1\tout2\tsingletons"
  echo -e "${SAMPLE}\t${pairs:-NA}\t${singletons:-NA}\t${time_s:-NA}\t${OUT1}\t${OUT2}\t${SING}"
} > "$STAT"

echo "[$(date)] Done ${SAMPLE}"
