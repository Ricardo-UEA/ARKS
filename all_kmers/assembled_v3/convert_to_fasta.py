# Convert kmers in plain text to FASTA format
with open("ARKS_master_kmers_filtered_2nd_screen.txt", "r") as input_file:
    with open("ARKS_master_kmers_filtered_2nd_screen.fasta", "w") as output_file:
        kmer_id = 1  # Start with a unique identifier for each kmer
        for line in input_file:
            line = line.strip()  # Remove any extra whitespace or newlines
            if line:  # Check if line is not empty
                # Write in FASTA format
                output_file.write(f">kmer{kmer_id}\n")  # Sequence ID
                output_file.write(f"{line}\n")          # Kmer sequence
                kmer_id += 1
