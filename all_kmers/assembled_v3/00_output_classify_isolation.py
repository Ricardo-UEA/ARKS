with open("ARKS_1st_assembly_kraken_output", "r") as infile:
    with open("unclassified_kmers_kraken.tsv", "w") as ufile, open("classified_kmers_kraken.tsv", "w") as cfile:
        for line in infile:
            if line.startswith("U"):
                ufile.write(line)
            elif line.startswith("C"):
                cfile.write(line)

print("Done: Saved unclassified_kmers.txt and classified_kmers.txt")
