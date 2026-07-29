#!/bin/bash
#SBATCH --job-name=camisim_simulations
#SBATCH --mail-type=ALL
#SBATCH --mail-user=hce24xau@uea.ac.uk
#SBATCH -o simulation_array_%j.out
#SBATCH -e simulation_array_%j.err
#SBATCH -p compute
#SBATCH --mem=400GB
#SBATCH --time=168:00:00
#SBATCH --cpus-per-task=20
#SBATCH --export=ALL
#SBATCH --array=1-1

module load apptainer

IMG=/gpfs/home/hce24xau/scratch/gen_kmers/data/simulations/CAMISIM/camisim.sif
WORKDIR=/gpfs/home/hce24xau/scratch/gen_kmers/data/simulations/CAMISIM

cd $WORKDIR

CONFIG=$WORKDIR/meta_simulations/configs/camisim_config_${SLURM_ARRAY_TASK_ID}.ini

apptainer exec --cleanenv --bind /gpfs:/gpfs camisim.sif 

python /gpfs/home/hce24xau/scratch/gen_kmers/data/simulations/CAMISIM/metagenomesimulation.py \
       $CONFIG \
       -p 20

echo "Processing meta_sample_${SLURM_ARRAY_TASK_ID}..."

SAMPLE_ROOT="${WORKDIR}/output_files/meta_sample_${SLURM_ARRAY_TASK_ID}"
SAMPLE_DIR=$(find "${SAMPLE_ROOT}" -maxdepth 1 -type d -name '*_sample_0' | head -n 1)

    # -------------------------
    # 1) Move gsa assembly
    # -------------------------
cd "${SAMPLE_DIR}/contigs"
mv gsa.fasta.gz "${WORKDIR}/meta_simulations/meta_assemblies/${SLURM_ARRAY_TASK_ID}.gsa.fasta.gz"

    # -------------------------
    # 2) Pool reads (R1 + R2)
    # -------------------------
cd "${SAMPLE_DIR}/reads"
zcat Genome*1.fq.gz | gzip > "${SAMPLE_DIR}/meta_sample_${SLURM_ARRAY_TASK_ID}_R1.fastq.gz"
zcat Genome*2.fq.gz | gzip > "${SAMPLE_DIR}/meta_sample_${SLURM_ARRAY_TASK_ID}_R2.fastq.gz"

#    # Move pooled reads
mv "${SAMPLE_DIR}/meta_sample_${SLURM_ARRAY_TASK_ID}_R1.fastq.gz" "${WORKDIR}/meta_simulations/meta_reads/"
mv "${SAMPLE_DIR}/meta_sample_${SLURM_ARRAY_TASK_ID}_R2.fastq.gz" "${WORKDIR}/meta_simulations/meta_reads/"

    # -------------------------
    # 3) Move taxonomic profile
    # -------------------------
mv "${SAMPLE_ROOT}/taxonomic_profile_0.txt" "${WORKDIR}/meta_simulations/taxonomic_profiles/meta_sample_${SLURM_ARRAY_TASK_ID}_taxonomic_profile.txt"

    # -------------------------
    # 4) Move distribution file
    # -------------------------
mv "${SAMPLE_ROOT}/distributions/distribution_0.txt" \
       "${WORKDIR}/meta_simulations/distributions/meta_sample_${SLURM_ARRAY_TASK_ID}_distribution.txt"

echo "All 1 ^`^s10 samples processed."
