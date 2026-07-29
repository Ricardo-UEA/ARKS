#!/bin/bash
#SBATCH --job-name=camisim_simulations
#SBATCH --mail-type=ALL
#SBATCH --mail-user=hce24xau@uea.ac.uk
#SBATCH -o job_outputs/simulation_array_%A_%a.out
#SBATCH -e job_errors/simulation_array_%A_%a.err
#SBATCH -p compute
#SBATCH --mem=400GB
#SBATCH --time=168:00:00
#SBATCH --cpus-per-task=1
#SBATCH --export=ALL
#SBATCH --array=1-100

module load apptainer
module load python/anaconda/2020.11

IMG=/gpfs/home/hce24xau/scratch/gen_kmers/data/simulations/CAMISIM/camisim.sif
WORKDIR=/gpfs/home/hce24xau/scratch/gen_kmers/data/simulations/CAMISIM

cd $WORKDIR

CONFIG=$WORKDIR/meta_simulations_2026_v2/configs/config_${SLURM_ARRAY_TASK_ID}.ini

echo "Starting meta_sample_${SLURM_ARRAY_TASK_ID}..."
echo "Using config: $CONFIG"

apptainer exec --cleanenv --bind /gpfs:/gpfs $IMG \
    python /gpfs/home/hce24xau/scratch/gen_kmers/data/simulations/CAMISIM/metagenomesimulation.py \
    $CONFIG \
    -p 1

echo "Finished meta_sample_${SLURM_ARRAY_TASK_ID}"
