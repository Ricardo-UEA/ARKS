#!/usr/bin/env python3
"""
Evaluate KrakenUniq classification on CAMISIM assemblies.

Compares assembly-based vs read-based classification.

Usage:
    python evaluate_assemblies.py /path/to/meta_simulations_2026 --samples 100 --output assembly_eval
"""

import argparse
import pandas as pd
import numpy as np
from pathlib import Path
from collections import defaultdict


def parse_taxonomic_profile(filepath):
    """Parse CAMISIM taxonomic_profile_0.txt and extract genus-level data."""
    genera = {}
    with open(filepath, 'r') as f:
        for line in f:
            if line.startswith('@') or line.startswith('#') or not line.strip():
                continue
            parts = line.strip().split('\t')
            if len(parts) >= 5:
                rank = parts[1]
                if rank == 'genus':
                    taxpath_sn = parts[3]
                    genus_name = taxpath_sn.split('|')[-1].strip()
                    percentage = float(parts[4])
                    if genus_name and percentage > 0:
                        genera[genus_name] = percentage
    return genera


def parse_kraken_report(filepath):
    """Parse KrakenUniq report for assemblies (contigs instead of reads)."""
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
                        taxname = parts[8].strip()
                        if taxname:
                            genera[taxname] = {
                                'pct': float(parts[0]),
                                'contigs': int(parts[1]),      # reads column = contigs for assembly
                                'taxContigs': int(parts[2]),   # taxReads = directly assigned contigs
                                'kmers': int(parts[3]),
                                'dup': float(parts[4]),
                                'cov': float(parts[5]),
                            }
                except (ValueError, IndexError):
                    continue
    return genera


def calculate_metrics_binary(truth_set, pred_set):
    """Calculate presence/absence metrics."""
    tp = len(truth_set & pred_set)
    fp = len(pred_set - truth_set)
    fn = len(truth_set - pred_set)
    
    precision = tp / (tp + fp) if (tp + fp) > 0 else 0
    recall = tp / (tp + fn) if (tp + fn) > 0 else 0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0
    
    return {'TP': tp, 'FP': fp, 'FN': fn, 'Precision': precision, 'Recall': recall, 'F1': f1}


def main():
    parser = argparse.ArgumentParser(description='Evaluate assembly-based classification')
    parser.add_argument('base_dir', help='Base directory')
    parser.add_argument('--samples', type=int, default=100, help='Number of samples')
    parser.add_argument('--output', default='assembly_evaluation', help='Output prefix')
    args = parser.parse_args()
    
    base_dir = Path(args.base_dir)
    assembly_reports_dir = base_dir / "kraken_assembly_reports"
    read_reports_dir = base_dir / "kraken_reports"
    output_files_dir = base_dir / "output_files"
    
    results = []
    
    print(f"Evaluating {args.samples} samples")
    print(f"Assembly reports: {assembly_reports_dir}")
    print(f"Read reports: {read_reports_dir}")
    print()
    
    for i in range(1, args.samples + 1):
        sample_name = f"meta_sample_{i}"
        
        # Files
        truth_file = output_files_dir / sample_name / "taxonomic_profile_0.txt"
        assembly_report = assembly_reports_dir / f"{sample_name}_assembly.report"
        read_report = read_reports_dir / f"{sample_name}.report"
        
        if not truth_file.exists():
            continue
        
        truth = parse_taxonomic_profile(truth_file)
        truth_set = set(truth.keys())
        
        row = {'Sample': sample_name, 'Truth_Genera': len(truth_set)}
        
        # Assembly-based classification
        if assembly_report.exists():
            assembly_pred = parse_kraken_report(assembly_report)
            
            # Test different contig thresholds (SEPATH used >=5)
            for min_contigs in [0, 1, 3, 5, 10]:
                pred_set = set(k for k, v in assembly_pred.items() if v['contigs'] >= min_contigs)
                m = calculate_metrics_binary(truth_set, pred_set)
                row[f'Assembly_Contigs>={min_contigs}_Prec'] = m['Precision']
                row[f'Assembly_Contigs>={min_contigs}_Rec'] = m['Recall']
                row[f'Assembly_Contigs>={min_contigs}_F1'] = m['F1']
        
        # Read-based classification (for comparison)
        if read_report.exists():
            read_pred = parse_kraken_report(read_report)
            
            # Best filter from previous analysis
            pred_set = set(k for k, v in read_pred.items() 
                          if v.get('kmers', 0) >= 100 and v.get('cov', 0) >= 0.002)
            m = calculate_metrics_binary(truth_set, pred_set)
            row['Reads_kmers100_cov002_Prec'] = m['Precision']
            row['Reads_kmers100_cov002_Rec'] = m['Recall']
            row['Reads_kmers100_cov002_F1'] = m['F1']
        
        results.append(row)
        
        if i % 10 == 0:
            print(f"  Processed {i}/{args.samples} samples")
    
    df = pd.DataFrame(results)
    
    # Summary
    print("\n" + "="*70)
    print("COMPARISON: Assembly-based vs Read-based Classification")
    print("="*70)
    
    print(f"\n{'Method':<40} {'Precision':>10} {'Recall':>10} {'F1':>10}")
    print("-"*70)
    
    # Assembly results
    for min_contigs in [0, 1, 3, 5, 10]:
        col = f'Assembly_Contigs>={min_contigs}_F1'
        if col in df.columns:
            prec = df[f'Assembly_Contigs>={min_contigs}_Prec'].mean()
            rec = df[f'Assembly_Contigs>={min_contigs}_Rec'].mean()
            f1 = df[col].mean()
            f1_std = df[col].std()
            label = f"Assembly (contigs >= {min_contigs})"
            print(f"{label:<40} {prec:>10.4f} {rec:>10.4f} {f1:>10.4f} +/- {f1_std:.4f}")
    
    # Read results
    if 'Reads_kmers100_cov002_F1' in df.columns:
        prec = df['Reads_kmers100_cov002_Prec'].mean()
        rec = df['Reads_kmers100_cov002_Rec'].mean()
        f1 = df['Reads_kmers100_cov002_F1'].mean()
        f1_std = df['Reads_kmers100_cov002_F1'].std()
        print(f"{'Reads (kmers>=100, cov>=0.002)':<40} {prec:>10.4f} {rec:>10.4f} {f1:>10.4f} +/- {f1_std:.4f}")
    
    # Save
    df.to_csv(f"{args.output}_results.csv", index=False)
    print(f"\nResults saved to {args.output}_results.csv")


if __name__ == '__main__':
    main()
