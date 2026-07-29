# Convert kmers in plain text to FASTA format
with open("T2T.txt", "r") as input_file:
    with open("T2T.fastq", "w") as output_file:
        kmer_id = 1  # Start with a unique identifier for each kmer
        for line in input_file:
            line = line.strip()  # Remove any extra whitespace or newlines
            if line:  # Check if line is not empty
                # Write in FASTA format
                output_file.write(f">kmer{kmer_id}\n")  # Sequence ID
                output_file.write(f"{line}\n")          # Kmer sequence
                kmer_id += 1
# Convert kmers in plain text to FASTA format
with open("GRCh37.txt", "r") as input_file:
    with open("GRCh37.fastq", "w") as output_file:
        kmer_id = 1  # Start with a unique identifier for each kmer
        for line in input_file:
            line = line.strip()  # Remove any extra whitespace or newlines
            if line:  # Check if line is not empty
                # Write in FASTA format
                output_file.write(f">kmer{kmer_id}\n")  # Sequence ID
                output_file.write(f"{line}\n")          # Kmer sequence
                kmer_id += 1
# Convert kmers in plain text to FASTA format
with open("GRCh38.txt", "r") as input_file:
    with open("GRCh38.fastq", "w") as output_file:
        kmer_id = 1  # Start with a unique identifier for each kmer
        for line in input_file:
            line = line.strip()  # Remove any extra whitespace or newlines
            if line:  # Check if line is not empty
                # Write in FASTA format
                output_file.write(f">kmer{kmer_id}\n")  # Sequence ID
                output_file.write(f"{line}\n")          # Kmer sequence
                kmer_id += 1

