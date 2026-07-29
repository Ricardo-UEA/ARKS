import duckdb

# File paths
csv_file = "/gpfs/home/hce24xau/scratch/gen_kmers/data/all_kmers/output_files/kmer_presence_table.csv"
db_file = "/gpfs/home/hce24xau/scratch/gen_kmers/data/all_kmers/output_files/kmer_presence.db"

# Connect to DuckDB (creates a new database if it doesn't exist)
con = duckdb.connect(db_file)

# Load CSV into DuckDB
con.execute(f"""
    CREATE OR REPLACE TABLE kmer_presence AS 
    SELECT * FROM read_csv_auto('{csv_file}', header=True);
""")

# Create an index on the 'kmer' column for faster lookups
con.execute("CREATE INDEX kmer_index ON kmer_presence (kmer);")

# Close the connection
con.close()

print(f"DuckDB database created successfully: {db_file}")
