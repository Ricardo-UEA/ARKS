import csv
import os

# Directory path for the input files
base_dir = "/gpfs/home/hce24xau/scratch/gen_kmers/data/"

# Define dataset paths
dataset_files = {
    "African_pan": base_dir + "african_pan/output_files/african_pan.txt",
    "GRCh37": base_dir + "references/output_files/GRCh37.txt",
    "GRCh38": base_dir + "references/output_files/GRCh38.txt",
    "T2T": base_dir + "references/output_files/T2T.txt",
    "HPRC2": base_dir + "pan_genome_2/output_files/pan_genome_2.txt",
    "IPD_IMGT_HLA": base_dir + "ipd_imgt_hla_db/output_files/hla_nuc_kmers.txt",
    "CPG": base_dir + "chinese_pan_genome/output_files/chinese_pan_genome.txt",
    "Arab-PG": base_dir + "arab_pan_genome/output_files/arab_pan_genome.txt",
    "KSPG": base_dir + "khoisan_pan_genome/output_files/KSPG.txt"
}

# Combined k-mer file
combined_kmers_file = base_dir + "all_kmers/assembled_v3/ARKS_master_kmers_filtered_2nd_screen.txt"

# Output file
output_csv_file = base_dir + "all_kmers/output_files/kmer_presence_table_all.csv"

# Function to load kmers from a file into a set
def load_kmers(file_path):
    with open(file_path, 'r') as file:
        return set(line.strip() for line in file)

# Load kmers from all dataset files
print("Loading kmers from dataset files...")
kmers_dict = {name: load_kmers(file) for name, file in dataset_files.items()}

# Generate the k-mer presence table
print("Generating k-mer presence table...")
with open(combined_kmers_file, 'r') as combined_file, open(output_csv_file, 'w', newline='') as output_file:
    csv_writer = csv.writer(output_file)
    
    # Write the header row
    csv_writer.writerow(["kmer"] + list(dataset_files.keys()))

    # Check presence of each k-mer in datasets
    for line in combined_file:
        kmer = line.strip()
        presence_flags = ["1" if kmer in kmers_dict[name] else "0" for name in dataset_files]
        csv_writer.writerow([kmer] + presence_flags)

print(f"K-mer presence table saved to {output_csv_file}")
