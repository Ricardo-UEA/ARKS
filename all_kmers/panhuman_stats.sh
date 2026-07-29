#!/bin/bash
#SBATCH --job-name=master_stats
#SBATCH -o master_stats_%j.out
#SBATCH -e master_stats_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=hce24xau@uea.ac.uk
#SBATCH -p hmem
#SBATCH --qos=hmem
#SBATCH --mem=500GB
#SBATCH --cpus-per-task=12
#SBATCH --time=168:00:00

module load python/anaconda/2020.11
source /gpfs/home/hce24xau/.bashrc
conda activate upset_plot

MASTER_FA="/gpfs/home/hce24xau/scratch/gen_kmers/data/all_kmers/assembled/ref_pans_hla_assembled.fasta"
OUTDIR="/gpfs/home/hce24xau/scratch/gen_kmers/data/all_kmers/results/master_stats"
K=31
JF_SIZE="40G"

mkdir -p "$OUTDIR"

MASTER_JF="${OUTDIR}/panHuman_master.jf"

jellyfish histo -t "$SLURM_CPUS_PER_TASK" -o "${OUTDIR}/kmer_histogram.txt" "$MASTER_JF"

awk '
BEGIN{total=0;distinct=0;max_freq=0;max_count=0}
{freq=$1;count=$2;total+=freq*count;distinct+=count;if(count>max_freq){max_freq=count;max_count=freq}}
END{
  print "Distinct_kmers\t" distinct;
  print "Total_kmer_obs\t" total;
  if(distinct>0) print "Mean_copy_number\t" total/distinct; else print "Mean_copy_number\tNA";
  print "Modal_copy_number\t" max_count;
}' "${OUTDIR}/kmer_histogram.txt" > "${OUTDIR}/kmer_histogram_summary.tsv"

awk '
{freq=$1;count=$2;
 if(freq==1) s+=count;
 else if(freq<=5) vr+=count;
 else if(freq<=50) m+=count;
 else if(freq<=500) c+=count;
 else hc+=count}
END{
 print "Singletons(=1)\t" (s+0);
 print "Very_rare(2-5)\t" (vr+0);
 print "Mid(6-50)\t" (m+0);
 print "Common(51-500)\t" (c+0);
 print "Highly_conserved(>500)\t" (hc+0);
}' "${OUTDIR}/kmer_histogram.txt" > "${OUTDIR}/kmer_copy_bins.tsv"

jellyfish dump -c -t "$MASTER_JF" | sort -k2,2rn | head -10 > "${OUTDIR}/top10_kmers.txt"

seqkit stats --all -N 10 -N 25 -N 50 -N 75 -N 90 --threads "$SLURM_CPUS_PER_TASK" "$MASTER_FA" > "${OUTDIR}/unitig_seqkit_stats.txt"

seqkit fx2tab --length --threads "$SLURM_CPUS_PER_TASK" "$MASTER_FA" | awk '{print $2}' | sort -n > "${OUTDIR}/unitig_lengths.txt"

awk '
{len=$1;total++;sum+=len;
 if(len==31) e31++;
 else if(len<=100) sh++;
 else if(len<=1000) me++;
 else lo++}
END{
 print "Total_unitigs\t" total;
 print "Total_bp\t" sum;
 if(total>0) print "Mean_len\t" sum/total; else print "Mean_len\tNA";
 print "Exact31\t" (e31+0);
 print "Short32-100\t" (sh+0);
 print "Medium101-1000\t" (me+0);
 print "Long>1000\t" (lo+0);
}' "${OUTDIR}/unitig_lengths.txt" > "${OUTDIR}/unitig_length_summary.tsv"

seqkit fx2tab --gc --threads "$SLURM_CPUS_PER_TASK" "$MASTER_FA" | awk '{print $3}' > "${OUTDIR}/gc_per_unitig.txt"

awk '
{sum+=$1;sumsq+=$1^2;n++;
 if($1<min||NR==1) min=$1;
 if($1>max) max=$1;
 bin[int($1)]++;
 if($1<35||$1>50) out++}
END{
 mean=sum/n; sd=sqrt(sumsq/n - mean^2);
 print "Mean_GC\t" mean;
 print "SD_GC\t" sd;
 print "Min_GC\t" min;
 print "Max_GC\t" max;
 print "Outliers(<35_or_>50)\t" (out+0);
 for(b=0;b<=100;b++) if(bin[b]>0) print "GC_" b "\t" bin[b];
}' "${OUTDIR}/gc_per_unitig.txt" > "${OUTDIR}/gc_summary.tsv"
