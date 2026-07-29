#!/bin/bash
#SBATCH --job-name=count_reads
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=hce24xau@uea.ac.uk
#SBATCH -o job_outputs/count_reads_%j.out
#SBATCH -e job_errors/count_reads_%j.err
#SBATCH -p hmem
#SBATCH --qos=hmem
#SBATCH --mem=360GB
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=16

BASE_DIR=/gpfs/home/hce24xau/scratch/gen_kmers/data/simulations/CAMISIM/meta_simulations_2026_v2
RAW_DIR=${BASE_DIR}/meta_reads
DEPLETED_DIR=${BASE_DIR}/depleted_reads
OUTPUT=${BASE_DIR}/depletion_read_counts.tsv

echo -e "sample\traw_reads\tdepleted_reads\tremoved_reads\tremoved_pct" > "$OUTPUT"

count_reads() {
    local file=$1
    echo $(zcat "$file" | wc -l) / 4 | bc
}

export -f count_reads

process_sample() {
    local sample=$1
    local RAW_DIR=$2
    local DEPLETED_DIR=$3

    R1_RAW="${RAW_DIR}/${sample}_R1.fq.gz"
    R1_DEP="${DEPLETED_DIR}/${sample}_depleted_R1.fq.gz"

    if [[ ! -f "$R1_RAW" ]]; then
        echo "${sample}: missing raw R1, skipping" >&2
        return
    fi
    if [[ ! -f "$R1_DEP" ]]; then
        echo "${sample}: missing depleted R1, skipping" >&2
        return
    fi

    RAW=$(zcat "$R1_RAW" | wc -l)
    RAW=$((RAW / 4))

    DEP=$(zcat "$R1_DEP" | wc -l)
    DEP=$((DEP / 4))

    REMOVED=$((RAW - DEP))
    PCT=$(echo "scale=4; $REMOVED * 100 / $RAW" | bc)

    echo -e "${sample}\t${RAW}\t${DEP}\t${REMOVED}\t${PCT}"
}

export -f process_sample

# Generate sample list
SAMPLES=$(for f in ${RAW_DIR}/meta_sample_*_R1.fq.gz; do
    basename "$f" _R1.fq.gz
done)

# Run in parallel
echo "$SAMPLES" | parallel -j 16 process_sample {} "$RAW_DIR" "$DEPLETED_DIR" >> "$OUTPUT"

# Sort by sample number
head -1 "$OUTPUT" > /tmp/header.tsv
tail -n +2 "$OUTPUT" | sort -t_ -k3 -n >> /tmp/sorted.tsv
cat /tmp/header.tsv /tmp/sorted.tsv > "$OUTPUT"

echo ""
echo "Done. Output written to: $OUTPUT"
echo ""
echo "Summary:"
awk -F'\t' 'NR>1 {sum_raw+=$2; sum_dep+=$3; sum_rem+=$4} 
     END {printf "Total raw reads:      %d\nTotal depleted reads: %d\nTotal removed:        %d\nOverall removal %%:    %.2f%%\n", 
     sum_raw, sum_dep, sum_rem, sum_rem*100/sum_raw}' "$OUTPUT"
