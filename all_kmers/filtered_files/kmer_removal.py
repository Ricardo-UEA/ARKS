import os

# Input microbial k-mer files (contaminants)
contaminant_files = [
    "/gpfs/home/hce24xau/scratch/gen_kmers/data/unique_kmers/african_pan/classified/microbial_kmers.txt",
    "/gpfs/home/hce24xau/scratch/gen_kmers/data/unique_kmers/pan_genome_2/classified/microbial_kmers.txt",
    "/gpfs/home/hce24xau/scratch/gen_kmers/data/unique_kmers/chinese_pan_genome/classified/microbial_kmers.txt"
]


# Master file (all kmers)
master_file = "/gpfs/home/hce24xau/scratch/gen_kmers/data/all_kmers/output_files/all_kmers_only_comb.txt"

# Output directory and file
output_dir = "/gpfs/home/hce24xau/scratch/gen_kmers/data/all_kmers/filtered_files"
os.makedirs(output_dir, exist_ok=True)
output_file = os.path.join(output_dir, "filtered_all_kmers.txt")

# Step 1: Combine all contaminant k-mers into a set
contaminant_kmers = set()
for filepath in contaminant_files:
    with open(filepath, "r") as f:
        contaminant_kmers.update(line.strip() for line in f if line.strip())

print(f"Total contaminant kmers loaded: {len(contaminant_kmers):,}")

# Step 2: Filter the master file
with open(master_file, "r") as infile, open(output_file, "w") as outfile:
    kept = 0
    for line in infile:
        kmer = line.strip()
        if kmer and kmer not in contaminant_kmers:
            outfile.write(kmer + "\n")
            kept += 1

print(f"Filtered output written to: {output_file}")
print(f"K-mers kept: {kept:,}")
