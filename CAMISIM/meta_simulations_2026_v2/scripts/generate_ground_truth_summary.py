#!/usr/bin/env python3
"""
Generate comprehensive ground truth summary for CAMISIM simulations.
Version 2: Uses NCBI taxonomy dump to map taxids to species names.

This script extracts:
- Read counts per genome per sample
- Taxonomic information (species, genus, family, etc.) from NCBI taxonomy
- Abundance percentages

Outputs:
- ground_truth_complete.csv: Complete summary of all samples
- ground_truth_by_sample.csv: Per-sample statistics
- ground_truth_by_taxon.csv: Aggregated statistics per taxon across all samples
- ground_truth_species_reads.csv: Simplified species view

Usage:
    python generate_ground_truth_summary_v2.py
"""

import os
import re
import gzip
import csv
from pathlib import Path
from collections import defaultdict
import subprocess

# ============================================================
# Configuration
# ============================================================

BASE_DIR = Path("/gpfs/home/hce24xau/scratch/gen_kmers/data/simulations/CAMISIM/meta_simulations_2026_v2")
OUTPUT_DIR = BASE_DIR / "output_files"
META_READS_DIR = BASE_DIR / "meta_reads"
METADATA_DIR = BASE_DIR / "metadata"
ID_TO_GENOME_DIR = BASE_DIR / "id_to_genomes"
MASTER_DIR = BASE_DIR / "master_files"
SUMMARY_DIR = BASE_DIR / "ground_truth_summary"

N_SAMPLES = 100

# ============================================================
# Hardcoded taxonomy mapping from master metadata
# This maps NCBI strain taxids to their taxonomic lineage
# ============================================================

# We'll build this from the taxonomic profile which has full lineages
TAXID_TO_TAXONOMY = {}

def build_taxonomy_from_profile(filepath: Path):
    """
    Parse a taxonomic profile to build taxid -> taxonomy mapping.
    We use strain-level entries which have complete lineage.
    """
    global TAXID_TO_TAXONOMY
    
    if not filepath.exists():
        return
    
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('@') or line.startswith('#') or not line:
                continue
            
            parts = line.split('\t')
            if len(parts) < 5:
                continue
            
            taxid = parts[0].split('.')[0]  # Remove .1 suffix
            rank = parts[1]
            taxpath = parts[2]
            taxpathsn = parts[3]
            
            # Parse the taxonomy path
            taxpath_parts = taxpathsn.split('|')
            
            # Store taxonomy for this taxid
            if taxid not in TAXID_TO_TAXONOMY:
                TAXID_TO_TAXONOMY[taxid] = {
                    'rank': rank,
                    'superkingdom': taxpath_parts[0] if len(taxpath_parts) > 0 else "",
                    'phylum': taxpath_parts[1] if len(taxpath_parts) > 1 else "",
                    'class': taxpath_parts[2] if len(taxpath_parts) > 2 else "",
                    'order': taxpath_parts[3] if len(taxpath_parts) > 3 else "",
                    'family': taxpath_parts[4] if len(taxpath_parts) > 4 else "",
                    'genus': taxpath_parts[5] if len(taxpath_parts) > 5 else "",
                    'species': taxpath_parts[6] if len(taxpath_parts) > 6 else "",
                    'strain': taxpath_parts[7] if len(taxpath_parts) > 7 else ""
                }


def count_reads_fastq(filepath: Path) -> int:
    """Count reads in a gzipped FASTQ file."""
    try:
        result = subprocess.run(
            f"zcat {filepath} | wc -l",
            shell=True,
            capture_output=True,
            text=True
        )
        lines = int(result.stdout.strip())
        return lines // 4
    except Exception as e:
        print(f"  Warning: Could not count reads in {filepath}: {e}")
        return 0


def parse_distribution_file(filepath: Path) -> dict:
    """Parse distribution_0.txt to get relative abundances per genome."""
    abundances = {}
    if not filepath.exists():
        return abundances
    
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) >= 2:
                genome_id = parts[0]
                abundance = float(parts[1])
                abundances[genome_id] = abundance
    
    return abundances


def parse_metadata_file(filepath: Path) -> dict:
    """Parse metadata TSV file to get NCBI IDs for each genome."""
    metadata = {}
    if not filepath.exists():
        return metadata
    
    with open(filepath, 'r') as f:
        header = f.readline().strip().split('\t')
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                genome_id = parts[0]
                otu = parts[1]
                ncbi_id = parts[2]
                novelty = parts[3] if len(parts) > 3 else ""
                metadata[genome_id] = {
                    'otu': otu,
                    'ncbi_id': ncbi_id,
                    'novelty_category': novelty
                }
    
    return metadata


def parse_id_to_genome_file(filepath: Path) -> dict:
    """Parse id_to_genome TSV file to get genome paths."""
    id_to_genome = {}
    if not filepath.exists():
        return id_to_genome
    
    with open(filepath, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                genome_id = parts[0]
                genome_path = parts[1]
                # Extract assembly accession from path
                assembly = Path(genome_path).stem.replace('_genomic', '')
                id_to_genome[genome_id] = {
                    'path': genome_path,
                    'assembly': assembly
                }
    
    return id_to_genome


def get_taxonomy_for_taxid(taxid: str) -> dict:
    """Get taxonomy info for a taxid from our built mapping."""
    taxid_clean = str(taxid).split('.')[0]
    
    if taxid_clean in TAXID_TO_TAXONOMY:
        return TAXID_TO_TAXONOMY[taxid_clean]
    
    return {
        'rank': '',
        'superkingdom': '',
        'phylum': '',
        'class': '',
        'order': '',
        'family': '',
        'genus': '',
        'species': '',
        'strain': ''
    }


def process_sample(sample_num: int) -> list:
    """Process a single sample and return list of genome records."""
    records = []
    
    dist_file = OUTPUT_DIR / f"meta_sample_{sample_num}" / "distributions" / "distribution_0.txt"
    tax_file = OUTPUT_DIR / f"meta_sample_{sample_num}" / "taxonomic_profile_0.txt"
    meta_file = METADATA_DIR / f"metadata_{sample_num}.tsv"
    id_genome_file = ID_TO_GENOME_DIR / f"id_to_genome_{sample_num}.tsv"
    reads_file = META_READS_DIR / f"meta_sample_{sample_num}_R1.fq.gz"
    
    # Build taxonomy from this sample's profile
    build_taxonomy_from_profile(tax_file)
    
    # Get total reads
    total_reads = 0
    if reads_file.exists():
        total_reads = count_reads_fastq(reads_file)
    
    # Parse files
    abundances = parse_distribution_file(dist_file)
    metadata = parse_metadata_file(meta_file)
    id_to_genome = parse_id_to_genome_file(id_genome_file)
    
    # Combine information for each genome
    for genome_id, abundance in abundances.items():
        reads = int(abundance * total_reads) if total_reads > 0 else 0
        abundance_pct = abundance * 100
        
        meta = metadata.get(genome_id, {})
        id_gen = id_to_genome.get(genome_id, {})
        
        # Get taxid from metadata
        taxid = meta.get('ncbi_id', '')
        
        # Get taxonomy from our mapping
        tax = get_taxonomy_for_taxid(taxid)
        
        record = {
            'sample': sample_num,
            'genome_id': genome_id,
            'assembly': id_gen.get('assembly', ''),
            'ncbi_taxid': taxid,
            'superkingdom': tax.get('superkingdom', ''),
            'phylum': tax.get('phylum', ''),
            'class': tax.get('class', ''),
            'order': tax.get('order', ''),
            'family': tax.get('family', ''),
            'genus': tax.get('genus', ''),
            'species': tax.get('species', ''),
            'strain': tax.get('strain', ''),
            'reads': reads,
            'abundance_pct': round(abundance_pct, 6),
            'abundance_fraction': abundance,
            'total_sample_reads': total_reads
        }
        records.append(record)
    
    return records


def main():
    print("=" * 60)
    print("CAMISIM Ground Truth Summary Generator v2")
    print("=" * 60)
    print(f"Base directory: {BASE_DIR}")
    print(f"Processing {N_SAMPLES} samples...")
    print()
    
    # Create output directory
    SUMMARY_DIR.mkdir(parents=True, exist_ok=True)
    
    # First pass: build taxonomy mapping from ALL samples
    print("Building taxonomy mapping from all samples...")
    for sample_num in range(1, N_SAMPLES + 1):
        tax_file = OUTPUT_DIR / f"meta_sample_{sample_num}" / "taxonomic_profile_0.txt"
        build_taxonomy_from_profile(tax_file)
    print(f"  Found {len(TAXID_TO_TAXONOMY)} unique taxids")
    print()
    
    all_records = []
    sample_stats = []
    taxon_stats = defaultdict(lambda: {
        'samples_present': 0,
        'total_reads': 0,
        'min_reads': float('inf'),
        'max_reads': 0,
        'min_abundance': float('inf'),
        'max_abundance': 0,
        'abundances': []
    })
    
    for sample_num in range(1, N_SAMPLES + 1):
        print(f"Processing sample {sample_num}/{N_SAMPLES}...", end=' ')
        
        records = process_sample(sample_num)
        all_records.extend(records)
        
        # Sample-level stats
        if records:
            total_reads = records[0]['total_sample_reads']
            n_genomes = len(records)
            
            sample_stats.append({
                'sample': sample_num,
                'total_reads': total_reads,
                'n_genomes': n_genomes,
                'min_abundance_pct': min(r['abundance_pct'] for r in records),
                'max_abundance_pct': max(r['abundance_pct'] for r in records),
                'min_reads': min(r['reads'] for r in records),
                'max_reads': max(r['reads'] for r in records)
            })
            
            # Taxon-level aggregation (by genus + species)
            for r in records:
                key = (r['genus'], r['species'])
                taxon_stats[key]['samples_present'] += 1
                taxon_stats[key]['total_reads'] += r['reads']
                taxon_stats[key]['min_reads'] = min(taxon_stats[key]['min_reads'], r['reads'])
                taxon_stats[key]['max_reads'] = max(taxon_stats[key]['max_reads'], r['reads'])
                taxon_stats[key]['min_abundance'] = min(taxon_stats[key]['min_abundance'], r['abundance_pct'])
                taxon_stats[key]['max_abundance'] = max(taxon_stats[key]['max_abundance'], r['abundance_pct'])
                taxon_stats[key]['abundances'].append(r['abundance_pct'])
                taxon_stats[key]['ncbi_taxid'] = r['ncbi_taxid']
                taxon_stats[key]['genus'] = r['genus']
                taxon_stats[key]['species'] = r['species']
                taxon_stats[key]['family'] = r['family']
                taxon_stats[key]['phylum'] = r['phylum']
                taxon_stats[key]['strain'] = r['strain']
            
            print(f"{n_genomes} genomes, {total_reads:,} reads")
        else:
            print("No data")
    
    # Write main summary CSV
    print()
    print("Writing output files...")
    
    main_csv = SUMMARY_DIR / "ground_truth_complete.csv"
    fieldnames = [
        'sample', 'genome_id', 'assembly', 'ncbi_taxid',
        'superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'strain',
        'reads', 'abundance_pct', 'abundance_fraction', 'total_sample_reads'
    ]
    
    with open(main_csv, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        # Sort by sample, then by reads descending
        sorted_records = sorted(all_records, key=lambda x: (x['sample'], -x['reads']))
        writer.writerows(sorted_records)
    
    print(f"  ✓ {main_csv}")
    
    # Write sample summary CSV
    sample_csv = SUMMARY_DIR / "ground_truth_by_sample.csv"
    sample_fieldnames = ['sample', 'total_reads', 'n_genomes', 'min_abundance_pct', 'max_abundance_pct', 'min_reads', 'max_reads']
    
    with open(sample_csv, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=sample_fieldnames)
        writer.writeheader()
        writer.writerows(sample_stats)
    
    print(f"  ✓ {sample_csv}")
    
    # Write taxon summary CSV
    taxon_csv = SUMMARY_DIR / "ground_truth_by_taxon.csv"
    taxon_records = []
    for key, stats in taxon_stats.items():
        avg_abundance = sum(stats['abundances']) / len(stats['abundances']) if stats['abundances'] else 0
        taxon_records.append({
            'genus': stats['genus'],
            'species': stats['species'],
            'strain': stats.get('strain', ''),
            'ncbi_taxid': stats['ncbi_taxid'],
            'family': stats['family'],
            'phylum': stats['phylum'],
            'samples_present': stats['samples_present'],
            'total_reads_all_samples': stats['total_reads'],
            'avg_abundance_pct': round(avg_abundance, 6),
            'min_abundance_pct': round(stats['min_abundance'], 6),
            'max_abundance_pct': round(stats['max_abundance'], 6),
            'min_reads': stats['min_reads'],
            'max_reads': stats['max_reads']
        })
    
    taxon_records = sorted(taxon_records, key=lambda x: -x['samples_present'])
    taxon_fieldnames = [
        'genus', 'species', 'strain', 'ncbi_taxid', 'family', 'phylum',
        'samples_present', 'total_reads_all_samples',
        'avg_abundance_pct', 'min_abundance_pct', 'max_abundance_pct',
        'min_reads', 'max_reads'
    ]
    
    with open(taxon_csv, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=taxon_fieldnames)
        writer.writeheader()
        writer.writerows(taxon_records)
    
    print(f"  ✓ {taxon_csv}")
    
    # Write a simple species-only summary for quick reference
    species_csv = SUMMARY_DIR / "ground_truth_species_reads.csv"
    species_records = []
    for r in all_records:
        species_records.append({
            'sample': r['sample'],
            'species': r['species'],
            'strain': r['strain'],
            'genus': r['genus'],
            'ncbi_taxid': r['ncbi_taxid'],
            'reads': r['reads'],
            'abundance_pct': r['abundance_pct']
        })
    
    species_records = sorted(species_records, key=lambda x: (x['sample'], -x['reads']))
    
    with open(species_csv, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['sample', 'species', 'strain', 'genus', 'ncbi_taxid', 'reads', 'abundance_pct'])
        writer.writeheader()
        writer.writerows(species_records)
    
    print(f"  ✓ {species_csv}")
    
    # Print summary statistics
    print()
    print("=" * 60)
    print("SUMMARY STATISTICS")
    print("=" * 60)
    
    total_samples = len(sample_stats)
    total_records = len(all_records)
    unique_species = len(set((r['genus'], r['species']) for r in all_records))
    total_reads_all = sum(s['total_reads'] for s in sample_stats)
    
    print(f"Total samples processed: {total_samples}")
    print(f"Total genome-sample combinations: {total_records}")
    print(f"Unique species: {unique_species}")
    print(f"Total reads across all samples: {total_reads_all:,}")
    
    if sample_stats:
        avg_reads = total_reads_all / total_samples
        avg_genomes = sum(s['n_genomes'] for s in sample_stats) / total_samples
        print(f"Average reads per sample: {avg_reads:,.0f}")
        print(f"Average genomes per sample: {avg_genomes:.1f}")
    
    # Abundance range
    all_abundances = [r['abundance_pct'] for r in all_records]
    if all_abundances:
        print(f"Abundance range: {min(all_abundances):.6f}% - {max(all_abundances):.4f}%")
    
    # Check for missing taxonomy
    missing_tax = sum(1 for r in all_records if not r['species'])
    if missing_tax > 0:
        print(f"WARNING: {missing_tax} records missing species names")
    
    print()
    print(f"Output directory: {SUMMARY_DIR}")
    print()
    print("Files created:")
    print(f"  - ground_truth_complete.csv     (all data)")
    print(f"  - ground_truth_by_sample.csv    (per-sample stats)")
    print(f"  - ground_truth_by_taxon.csv     (per-taxon aggregated)")
    print(f"  - ground_truth_species_reads.csv (simplified species view)")


if __name__ == "__main__":
    main()
