#!/bin/bash
#SBATCH --job-name=process_kmers         # Job name
#SBATCH --mail-type=ALL                  # Mail notifications (ALL, BEGIN, END, FAIL)
#SBATCH --mail-user=hce24xau@uea.ac.uk   # Email address for notifications
#SBATCH -o process_kmers-%j.out          # Standard output log (%j for job ID)
#SBATCH -e process_kmers-%j.err          # Standard error log (%j for job ID)
#SBATCH -p hmem                      # Which queue to use
#SBATCH --qos=hmem                       # Access to the high-memory queue
#SBATCH --mem=2000GB                      # Memory required (adjust as needed)
#SBATCH --time=168:00:00                  # Time limit (adjust based on your job's estimated runtime)
#SBATCH --export=ALL                     # Export all environment variables
#SBATCH --ntasks=1                        # Run a single task (script)
#SBATCH --cpus-per-task=24   


# Load the required Python environment (or activate virtual environment)
module load python/anaconda/2020.11

source activate upset_plot

# Record the start time of the script
echo "Job started at $(date)"

# Change to the directory containing your Python scripts
cd /gpfs/home/hce24xau/scratch/gen_kmers/data/all_kmers || { 
  echo "Directory does not exist. Exiting."; exit 1; 
}

# Function to monitor memory usage
log_memory_usage() {
  echo "Logging memory usage to ${SLURM_JOB_ID}_memory.log"
  while true; do
    # Log memory usage with timestamp
    echo "$(date):" >> "${SLURM_JOB_ID}_memory.log"
    free -h >> "${SLURM_JOB_ID}_memory.log"
    echo "----------------------------------------" >> "${SLURM_JOB_ID}_memory.log"
    sleep 60  # Log memory usage every 60 seconds
  done
}

# Start memory monitoring in the background
log_memory_usage &
MEMORY_LOGGER_PID=$!  # Capture the PID of the memory logger

# Step 1: Process GRCh37 k-mers
echo "Processing GRCh37 k-mers..."
python generate_kmer_table.py

# Check the exit status of the Python script
if [ $? -eq 0 ]; then
  echo "Python script completed successfully."
else
  echo "Python script failed. Check error log for details." >&2
  kill $MEMORY_LOGGER_PID  # Stop the memory logger
  exit 1
fi

# Stop memory monitoring
kill $MEMORY_LOGGER_PID
wait $MEMORY_LOGGER_PID 2>/dev/null

# Record the end time of the script
echo "Job completed successfully at $(date)."
