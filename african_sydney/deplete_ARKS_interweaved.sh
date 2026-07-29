#!/bin/bash
#SBATCH --job-name=bbduk_african_reads_panhuman
#SBATCH --mail-type=ALL
#SBATCH --mail-user=hce24xau@uea.ac.uk
#SBATCH -o /gpfs/home/hce24xau/scratch/gen_kmers/data/african_sydney/job_outputs/bbduk_depletion/hdist_0_ARKS/bbduk_%A_%a.out
#SBATCH -e /gpfs/home/hce24xau/scratch/gen_kmers/data/african_sydney/job_errors/bbduk_depletion/hdist_0_ARKS/bbduk_%A_%a.err
#SBATCH -p compute
#SBATCH --mem=200GB
#SBATCH --time=24:00:00
#SBATCH --export=ALL
#SBATCH --cpus-per-task=6
#SBATCH --array=1-176

module load bbmap
module load python/anaconda/2020.11

BASE="/gpfs/home/hce24xau/scratch/gen_kmers/data/african_sydney"
READS_DIR="${BASE}/input_files/african_batch1"
SAMPLE_LIST="${READS_DIR}/samples.txt"

FILE_NAME="hdist_0_ARK_mcf_05"

OUT_BASE="${BASE}/output_files/${FILE_NAME}"
CLEAN_DIR="${OUT_BASE}/clean"
HUMAN_DIR="${OUT_BASE}/human"
STATS_DIR="${BASE}/stats_files/${FILE_NAME}"

#HUMAN_REF="/gpfs/home/hce24xau/scratch/gen_kmers/data/references/input_files/T2T.fasta"
#HUMAN_REF="/gpfs/home/hce24xau/scratch/gen_kmers/data/references/input_files/GRCh37.fasta"
#HUMAN_REF="/gpfs/home/hce24xau/scratch/gen_kmers/data/references/input_files/GRCh38.fasta"
HUMAN_REF="/gpfs/afm/hpccancergenetics/Ricardo/ARKS_assembly/ARKS_assembly.fasta"

mkdir -p "${CLEAN_DIR}" "${HUMAN_DIR}" "${STATS_DIR}" \
         "${BASE}/job_outputs" "${BASE}/job_errors/bbduk_depletion"

SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLE_LIST")

if [[ -z "$SAMPLE" ]]; then
    echo "No sample found for task ${SLURM_ARRAY_TASK_ID}" >&2
    exit 1
fi

R1="${READS_DIR}/${SAMPLE}_R1.fastq.gz"
R2="${READS_DIR}/${SAMPLE}_R2.fastq.gz"

echo "[$(date)] Starting: ${SAMPLE}"

if [[ ! -f "$R1" || ! -f "$R2" ]]; then
    echo "FASTQs missing for ${SAMPLE}" >&2
    echo "Expected: $R1 and $R2" >&2
    exit 1
fi

if [[ -f "${CLEAN_DIR}/${SAMPLE}_clean_R1.fastq.gz" && -f "${CLEAN_DIR}/${SAMPLE}_clean_R2.fastq.gz" ]]; then
    echo "Outputs already exist for ${SAMPLE}, skipping."
    exit 0
fi

bbduk.sh \
    -Xmx200g \
    in1="$R1" in2="$R2" \
    out="${CLEAN_DIR}/${SAMPLE}_clean.fastq.gz" \
    outm="${HUMAN_DIR}/${SAMPLE}_human.fastq.gz" \
    ref="$HUMAN_REF" \
    interleaved=t \
    k=31 \
    mcf=0.5 \
    hdist=0 \
    mm=t \
    removeifeitherbad=t \
    ordered=t \
    threads="$SLURM_CPUS_PER_TASK" \
    stats="${STATS_DIR}/${SAMPLE}_bbduk_stats.txt"


echo "[$(date)] Finished: ${SAMPLE}"
