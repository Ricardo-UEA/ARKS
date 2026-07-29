#!/bin/bash
#SBATCH --job-name=process_kmers         # Job name
#SBATCH --mail-type=ALL                  # Mail notifications (ALL, BEGIN, END, FAIL)
#SBATCH --mail-user=hce24xau@uea.ac.uk   # Email address for notifications
#SBATCH -o process_kmers-%j.out          # Standard output log (%j for job ID)
#SBATCH -e process_kmers-%j.err          # Standard error log (%j for job ID)
#SBATCH -p hmem                      # Which queue to use
#SBATCH --qos=hmem
#SBATCH --mem=2000GB                      # Memory required (adjust as needed)
#SBATCH --time=168:00:00                  # Time limit (adjust based on your job's estimated runtime)
#SBATCH --export=ALL                     # Export all environment variables
#SBATCH --ntasks=1                        # Run a single task (script)
#SBATCH --cpus-per-task=24   


module load python/anaconda/2020.11

source /gpfs/software/hali/python/anaconda/2020.11/3.8/etc/profile.d/conda.sh
conda activate upset_plot

echo "Python being used:"
which python

echo "Conda env:"
echo "$CONDA_DEFAULT_ENV"

pip install duckdb


python -c "import duckdb; print(duckdb.__version__)"

cd /gpfs/home/hce24xau/scratch/gen_kmers/data/all_kmers
python convert_to_duckdb.py
