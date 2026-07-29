#!/usr/bin/env python3
"""
Build genome x sample read count matrix - PARALLEL VERSION
"""

import gzip
from pathlib import Path
from collections import defaultdict
import pandas as pd
import re
from multiprocessing import Pool, cpu_count

BASE_DIR = Path("/gpfs/home/hce24xau/scratch/gen_kmers/data/simulations/CAMISIM/meta_simulations_2026_v2")
N_SAMPLES = 100
N_WORKERS = 16  # Adjust based on available cores


def find_reads_mapping(sample_num):
    sample_dir = BASE_DIR / f"output_files/meta_sample_{sample_num}"
    if not sample_dir.exists():
        return None
    for d in sample_dir.iterdir():
        if d.is_dir() and re.match(r'\d{4}\.\d{2}\.\d{2}', d.name):
            mapping_file = d / "reads" / "reads_mapping.tsv.gz"
            if mapping_file.exists():
                return mapping_file
    return None


def load_contig_to_genome_mapping(mapping_file):
    contig_map = {}
    with gzip.open(mapping_file, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                contig_map[parts[0]] = {'genome_id': parts[1], 'tax_id': parts[2]}
    return contig_map


def get_genome_to_assembly(sample_num):
    id_file = BASE_DIR / f"id_to_genomes/id_to_genome_{sample_num}.tsv"
    mapping = {}
    if not id_file.exists():
        return mapping
    with open(id_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                mapping[parts[0]] = Path(parts[1]).stem.replace('_genomic', '')
    return mapping


def process_sample(sample_num):
    """Process a single sample - designed for parallel execution."""
    mapping_file = find_reads_mapping(sample_num)
    if not mapping_file:
        return None, f"Sample {sample_num}: no mapping file"
    
    fastq_path = BASE_DIR / f"meta_reads/meta_sample_{sample_num}_R1.fq.gz"
    if not fastq_path.exists():
        return None, f"Sample {sample_num}: no FASTQ"
    
    contig_map = load_contig_to_genome_mapping(mapping_file)
    genome_to_assembly = get_genome_to_assembly(sample_num)
    
    # Get tax_id mapping
    genome_to_taxid = {v['genome_id']: v['tax_id'] for v in contig_map.values()}
    
    # Count reads
    genome_counts = defaultdict(int)
    total = 0
    unmapped = 0
    
    with gzip.open(fastq_path, 'rt') as f:
        for i, line in enumerate(f):
            if i % 4 == 0:
                total += 1
                contig_id = line.strip()[1:].rsplit('-', 1)[0]
                if contig_id in contig_map:
                    genome_counts[contig_map[contig_id]['genome_id']] += 1
                else:
                    unmapped += 1
    
    # Build results
    results = []
    for genome_id, read_count in genome_counts.items():
        results.append({
            'sample': sample_num,
            'genome_id': genome_id,
            'assembly': genome_to_assembly.get(genome_id, ''),
            'taxid': genome_to_taxid.get(genome_id, ''),
            'reads': read_count
        })
    
    msg = f"Sample {sample_num}: {len(genome_counts)} genomes, {total:,} reads, {unmapped:,} unmapped"
    return results, msg


def main():
    print("=" * 70)
    print(f"Building matrix from FASTQ headers - PARALLEL ({N_WORKERS} workers)")
    print("=" * 70)
    
    # Process in parallel
    with Pool(N_WORKERS) as pool:
        results = pool.map(process_sample, range(1, N_SAMPLES + 1))
    
    # Combine results
    all_data = []
    for data, msg in results:
        print(f"  {msg}")
        if data:
            all_data.extend(data)
    
    df = pd.DataFrame(all_data)
    
    # Save outputs
    output_dir = BASE_DIR / "ground_truth_summary"
    output_dir.mkdir(exist_ok=True)
    
    df.to_csv(output_dir / "reads_per_genome_from_fastq.csv", index=False)
    
    for idx_col, name in [('genome_id', 'genome'), ('taxid', 'taxid'), ('assembly', 'assembly')]:
        matrix = df.pivot_table(index=idx_col, columns='sample', values='reads', aggfunc='sum', fill_value=0)
        matrix = matrix[sorted(matrix.columns)]
        matrix.to_csv(output_dir / f"{name}_sample_reads_matrix.csv")
        print(f"Saved: {name}_sample_reads_matrix.csv ({matrix.shape})")
    
    print(f"\nTotal: {len(df)} records, {df['sample'].nunique()} samples, {df['taxid'].nunique()} taxids")


if __name__ == "__main__":
    main()
