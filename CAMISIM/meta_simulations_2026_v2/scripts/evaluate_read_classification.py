#!/usr/bin/env python3
"""
Evaluate KrakenUniq read-level classification accuracy against CAMISIM ground truth.

Metrics:
- Per-genus recall: classified_reads / true_reads
- Overall sensitivity: total_correctly_classified / total_true_reads
- False positive reads: reads assigned to genera not in ground truth

Usage:
    python evaluate_read_classification.py
"""

import pandas as pd
import numpy as np
from pathlib import Path
from collections import defaultdict

BASE_DIR = Path("/gpfs/home/hce24xau/scratch/gen_kmers/data/simulations/CAMISIM/meta_simulations_2026_v2")
N_SAMPLES = 100

# Taxid to genus mapping (same as before)
TAXID_TO_GENUS = {
    479436: "Veillonella",
    883161: "Vaginimicrobium",
    1123001: "Vaginimicrobium",
    1111131: "Vaginimicrobium",
    565575: "Ureaplasma",
    1125702: "Treponema",
    243276: "Treponema",
    428126: "Thomasclavelia",
    365659: "Streptococcus",
    760570: "Streptococcus",
    1051074: "Streptococcus",
    862971: "Streptococcus",
    1311: "Streptococcus",
    1309: "Streptococcus",
    519441: "Streptobacillus",
    523796: "Staphylococcus",
    1282: "Staphylococcus",
    40543: "Sneathia",
    10298: "Simplexvirus",
    546271: "Selenomonas",
    883077: "Schaalia",
    220341: "Salmonella",
    762948: "Rothia",
    37296: "Rhadinovirus",
    329: "Ralstonia",
    714315: "Pseudoleptotrichia",
    868129: "Prevotella",
    1122981: "Prevotella",
    1235811: "Prevotella",
    431947: "Porphyromonas",
    879243: "Porphyromonas",
    1122975: "Porphyromonas",
    1122971: "Porphyromonas",
    887901: "Porphyromonas",
    1583331: "Porphyromonas",
    28124: "Porphyromonas",
    2811780: "Porphyromonas",
    5858: "Plasmodium",
    862517: "Peptoniphilus",
    3036303: "Peptoniphilus",
    2811779: "Peptoniphilus",
    667128: "Pasteurella",
    11082: "Orthoflavivirus",
    435832: "Neisseria",
    546263: "Neisseria",
    495: "Neisseria",
    2097: "Mycoplasmoides",
    83332: "Mycobacterium",
    1236608: "Moraxella",
    114527: "Mogibacterium",
    548479: "Mobiluncus",
    411460: "Mediterraneibacter",
    1316932: "Mannheimia",
    10376: "Lymphocryptovirus",
    169963: "Listeria",
    47671: "Lautropia",
    553184: "Lancefieldella",
    629741: "Kingella",
    679190: "Hoylesella",
    1122992: "Hoylesella",
    41856: "Hepacivirus",
    85962: "Helicobacter",
    71421: "Haemophilus",
    46124: "Granulicatella",
    546270: "Gemella",
    356663: "Gammaretrovirus",
    190304: "Fusobacterium",
    469605: "Fusobacterium",
    525282: "Finegoldia",
    411483: "Faecalibacterium",
    411463: "Eubacterium",
    511145: "Escherichia",
    226185: "Enterococcus",
    1169293: "Enterococcus",
    164759: "Diaphorobacter",
    592028: "Dialister",
    1654930: "Cytomegalovirus",
    267747: "Cutibacterium",
    12305: "Cucumovirus",
    548477: "Corynebacterium",
    525264: "Corynebacterium",
    525260: "Corynebacterium",
    553206: "Corynebacterium",
    1125779: "Corynebacterium",
    61592: "Corynebacterium",
    45242: "Capnocytophaga",
    553218: "Campylobacter",
    1121102: "Campylobacter",
    1032069: "Campylobacter",
    199: "Campylobacter",
    679192: "Bulleidia",
    537007: "Blautia",
    1891762: "Betapolyomavirus",
    10632: "Betapolyomavirus",
    295405: "Bacteroides",
    326423: "Bacillus",
    879305: "Anaerococcus",
    525919: "Anaerococcus",
    561177: "Anaerococcus",
    1284686: "Anaerococcus",
    626522: "Alloprevotella",
    1120933: "Actinotignum",
    59505: "Actinotignum",
    435830: "Actinomyces",
    544580: "Actinomyces",
    592010: "Abiotrophia",
}

GROUND_TRUTH_GENERA = set(TAXID_TO_GENUS.values())


def parse_kraken_report(filepath):
    """
    Parse KrakenUniq report and extract genus-level read counts.
    Returns dict: genus_name -> reads
    """
    genera = {}
    
    with open(filepath, 'r') as f:
        for line in f:
            if line.startswith('#') or line.startswith('%'):
                continue
            parts = line.strip().split('\t')
            if len(parts) >= 9:
                try:
                    rank = parts[7]
                    if rank == 'genus':
                        reads = int(parts[1])      # Total reads at this clade
                        tax_reads = int(parts[2])  # Reads directly assigned
                        genus_name = parts[8].strip()
                        kmers = int(parts[3])
                        cov = float(parts[5])
                        
                        genera[genus_name] = {
                            'reads': reads,
                            'taxReads': tax_reads,
                            'kmers': kmers,
                            'cov': cov
                        }
                except (ValueError, IndexError):
                    continue
    
    return genera


def load_ground_truth():
    """Load genus-level ground truth from CSV."""
    gt_file = BASE_DIR / "ground_truth_summary" / "genus_sample_reads_matrix.csv"
    
    if not gt_file.exists():
        # Try to build from reads_per_genome_from_fastq.csv
        raw_file = BASE_DIR / "ground_truth_summary" / "reads_per_genome_from_fastq.csv"
        df = pd.read_csv(raw_file)
        df['genus'] = df['taxid'].map(TAXID_TO_GENUS)
        
        gt = df.pivot_table(
            index='genus',
            columns='sample',
            values='reads',
            aggfunc='sum',
            fill_value=0
        )
    else:
        gt = pd.read_csv(gt_file, index_col=0)
    
    return gt


def main():
    print("=" * 70)
    print("KrakenUniq Read-Level Classification Evaluation")
    print("=" * 70)
    
    # Load ground truth
    print("\nLoading ground truth...")
    gt_matrix = load_ground_truth()
    print(f"  Ground truth: {gt_matrix.shape[0]} genera x {gt_matrix.shape[1]} samples")
    
    # Store results
    all_results = []
    per_genus_stats = defaultdict(lambda: {'true_reads': 0, 'classified_reads': 0, 'samples': 0})
    
    print(f"\nProcessing {N_SAMPLES} Kraken reports...")
    
    for sample_num in range(1, N_SAMPLES + 1):
        report_file = BASE_DIR / "kraken_reports" / f"meta_sample_{sample_num}.report"
        
        if not report_file.exists():
            print(f"  Sample {sample_num}: report not found")
            continue
        
        # Get ground truth for this sample
        sample_col = str(sample_num) if str(sample_num) in gt_matrix.columns else sample_num
        if sample_col not in gt_matrix.columns:
            continue
            
        gt_sample = gt_matrix[sample_col]
        true_genera = set(gt_sample[gt_sample > 0].index)
        
        # Parse Kraken report
        kraken_genera = parse_kraken_report(report_file)
        
        # Calculate metrics for this sample
        sample_true_reads = 0
        sample_classified_reads = 0
        sample_fp_reads = 0
        
        # For each true genus, check if classified
        for genus in true_genera:
            true_reads = gt_sample[genus]
            sample_true_reads += true_reads
            
            if genus in kraken_genera:
                classified = kraken_genera[genus]['reads']
                sample_classified_reads += min(classified, true_reads)  # Cap at true reads
                
                per_genus_stats[genus]['true_reads'] += true_reads
                per_genus_stats[genus]['classified_reads'] += classified
                per_genus_stats[genus]['samples'] += 1
            else:
                per_genus_stats[genus]['true_reads'] += true_reads
                per_genus_stats[genus]['samples'] += 1
        
        # Count FP reads (assigned to genera not in ground truth)
        for genus, data in kraken_genera.items():
            if genus not in true_genera:
                sample_fp_reads += data['reads']
        
        sample_recall = sample_classified_reads / sample_true_reads if sample_true_reads > 0 else 0
        
        all_results.append({
            'sample': sample_num,
            'true_reads': sample_true_reads,
            'classified_reads': sample_classified_reads,
            'fp_reads': sample_fp_reads,
            'recall': sample_recall,
            'true_genera': len(true_genera),
            'detected_genera': len([g for g in true_genera if g in kraken_genera])
        })
        
        if sample_num % 20 == 0:
            print(f"  Processed {sample_num}/{N_SAMPLES} samples")
    
    # Create results DataFrame
    results_df = pd.DataFrame(all_results)
    
    # Summary statistics
    print("\n" + "=" * 70)
    print("OVERALL RESULTS")
    print("=" * 70)
    
    print(f"\nSamples processed: {len(results_df)}")
    print(f"\nRead-level metrics (mean ± std across samples):")
    print(f"  True reads per sample:       {results_df['true_reads'].mean():,.0f} ± {results_df['true_reads'].std():,.0f}")
    print(f"  Classified reads per sample: {results_df['classified_reads'].mean():,.0f} ± {results_df['classified_reads'].std():,.0f}")
    print(f"  FP reads per sample:         {results_df['fp_reads'].mean():,.0f} ± {results_df['fp_reads'].std():,.0f}")
    print(f"  Read recall:                 {results_df['recall'].mean():.4f} ± {results_df['recall'].std():.4f}")
    
    print(f"\nGenus-level detection:")
    print(f"  True genera per sample:      {results_df['true_genera'].mean():.1f}")
    print(f"  Detected genera per sample:  {results_df['detected_genera'].mean():.1f}")
    print(f"  Detection rate:              {results_df['detected_genera'].mean() / results_df['true_genera'].mean():.4f}")
    
    # Per-genus breakdown
    print("\n" + "=" * 70)
    print("PER-GENUS READ RECALL")
    print("=" * 70)
    
    genus_results = []
    for genus, stats in per_genus_stats.items():
        recall = stats['classified_reads'] / stats['true_reads'] if stats['true_reads'] > 0 else 0
        genus_results.append({
            'genus': genus,
            'total_true_reads': stats['true_reads'],
            'total_classified_reads': stats['classified_reads'],
            'recall': recall,
            'samples_present': stats['samples']
        })
    
    genus_df = pd.DataFrame(genus_results).sort_values('total_true_reads', ascending=False)
    
    print(f"\n{'Genus':<25} {'True Reads':>15} {'Classified':>15} {'Recall':>10} {'Samples':>10}")
    print("-" * 80)
    for _, row in genus_df.iterrows():
        print(f"{row['genus']:<25} {row['total_true_reads']:>15,} {row['total_classified_reads']:>15,} {row['recall']:>10.4f} {row['samples_present']:>10}")
    
    # Save results
    output_dir = BASE_DIR / "ground_truth_summary"
    results_df.to_csv(output_dir / "read_classification_per_sample.csv", index=False)
    genus_df.to_csv(output_dir / "read_classification_per_genus.csv", index=False)
    
    print(f"\nResults saved to:")
    print(f"  - read_classification_per_sample.csv")
    print(f"  - read_classification_per_genus.csv")


if __name__ == "__main__":
    main()
