#!/usr/bin/env python3
"""
Evaluate metagenomic classification performance across CAMISIM simulations.

Optimised filters based on grid search:
  - min_kmers >= 100
  - min_cov >= 0.001

Usage:
    python evaluate_simulations.py /gpfs/home/hce24xau/scratch/gen_kmers/data/simulations/CAMISIM/meta_simulations_2026_v2 --samples 100 --output evaluate_simulations_results/results

HPC Directory structure:
    base_dir/
    ├── kraken_reports/meta_sample_1.report ... meta_sample_N.report
    └── output_files/meta_sample_1/taxonomic_profile_0.txt ... meta_sample_N/taxonomic_profile_0.txt

Local (all in one directory):
    base_dir/
    ├── taxonomic_profile_1.txt ... taxonomic_profile_N.txt
    └── meta_sample_1.report ... meta_sample_N.report
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
    """Parse KrakenUniq report and extract ALL metrics for genus-level."""
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
                                'reads': int(parts[1]),
                                'taxReads': int(parts[2]),
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
    
    return {
        'TP': tp, 'FP': fp, 'FN': fn,
        'Precision_PPV': precision,
        'Recall_Sensitivity': recall,
        'F1': f1
    }


def calculate_metrics_abundance(truth_dict, pred_dict):
    """Calculate abundance-weighted metrics."""
    all_genera = set(truth_dict.keys()) | set(pred_dict.keys())
    
    if not all_genera:
        return {'Bray_Curtis_Dissimilarity': 1.0, 'Abundance_Correlation': 0.0, 'L1_Error': 0.0}
    
    truth_total = sum(truth_dict.values()) if truth_dict else 1
    pred_total = sum(pred_dict.values()) if pred_dict else 1
    
    truth_vec = [truth_dict.get(g, 0) / truth_total for g in all_genera]
    pred_vec = [pred_dict.get(g, 0) / pred_total for g in all_genera]
    
    truth_arr = np.array(truth_vec)
    pred_arr = np.array(pred_vec)
    
    bc = np.sum(np.abs(truth_arr - pred_arr)) / (np.sum(truth_arr) + np.sum(pred_arr)) if (np.sum(truth_arr) + np.sum(pred_arr)) > 0 else 0
    
    if np.std(truth_arr) > 0 and np.std(pred_arr) > 0:
        corr = np.corrcoef(truth_arr, pred_arr)[0, 1]
    else:
        corr = 0
    
    l1 = np.mean(np.abs(truth_arr - pred_arr))
    
    return {'Bray_Curtis_Dissimilarity': bc, 'Abundance_Correlation': corr, 'L1_Error': l1}


def find_files(base_dir, sample_num):
    """Find truth and report files for a sample."""
    base_dir = Path(base_dir)
    sample_name = f"meta_sample_{sample_num}"
    
    # Truth file locations
    truth_candidates = [
        base_dir / f"taxonomic_profile_{sample_num}.txt",
        base_dir / f"taxonomic_profile_{sample_num-1}.txt",
        base_dir / sample_name / "taxonomic_profile_0.txt",
        base_dir / "output_files" / sample_name / "taxonomic_profile_0.txt",
    ]
    
    # Report file locations
    report_candidates = [
        base_dir / f"{sample_name}.report",
        base_dir / f"meta_sample_{sample_num}.report",
        base_dir / "kraken_reports" / f"{sample_name}.report",
    ]
    
    truth_file = next((f for f in truth_candidates if f.exists()), None)
    report_file = next((f for f in report_candidates if f.exists()), None)
    
    return truth_file, report_file


def apply_filter(pred_data, min_reads=0, min_kmers=0, min_cov=0, min_taxReads=0, max_dup=float('inf')):
    """Apply filtering criteria to predictions."""
    return {
        k: v for k, v in pred_data.items()
        if v['reads'] >= min_reads
        and v['kmers'] >= min_kmers
        and v['cov'] >= min_cov
        and v['taxReads'] >= min_taxReads
        and v['dup'] <= max_dup
    }


def test_threshold(all_pred_full, truth_matrix, samples, **filter_kwargs):
    """Test a specific threshold combination across all samples."""
    metrics_list = []
    for i in range(1, samples + 1):
        sample_name = f"meta_sample_{i}"
        if sample_name not in truth_matrix.index or sample_name not in all_pred_full:
            continue
        truth_set = set(truth_matrix.columns[truth_matrix.loc[sample_name] > 0])
        pred_filtered = apply_filter(all_pred_full[sample_name], **filter_kwargs)
        pred_set = set(pred_filtered.keys())
        m = calculate_metrics_binary(truth_set, pred_set)
        metrics_list.append(m)
    
    if not metrics_list:
        return None
    
    return {
        'Precision': np.mean([m['Precision_PPV'] for m in metrics_list]),
        'Precision_std': np.std([m['Precision_PPV'] for m in metrics_list]),
        'Recall': np.mean([m['Recall_Sensitivity'] for m in metrics_list]),
        'Recall_std': np.std([m['Recall_Sensitivity'] for m in metrics_list]),
        'F1': np.mean([m['F1'] for m in metrics_list]),
        'F1_std': np.std([m['F1'] for m in metrics_list]),
    }


def main():
    parser = argparse.ArgumentParser(description='Evaluate CAMISIM simulation classification')
    parser.add_argument('base_dir', help='Directory containing files')
    parser.add_argument('--samples', type=int, default=100, help='Number of samples (default: 100)')
    parser.add_argument('--output', default='simulation_evaluation', help='Output prefix')
    # Optimised defaults from grid search
    parser.add_argument('--min-kmers', type=int, default=100, help='Minimum k-mers (default: 100)')
    parser.add_argument('--min-cov', type=float, default=0.001, help='Minimum coverage (default: 0.001)')
    parser.add_argument('--min-reads', type=int, default=0, help='Minimum reads (default: 0)')
    parser.add_argument('--min-taxReads', type=int, default=0, help='Minimum taxReads (default: 0)')
    parser.add_argument('--max-dup', type=float, default=100, help='Maximum duplication (default: 100)')
    args = parser.parse_args()
    
    base_dir = Path(args.base_dir)
    
    all_truth_genera = defaultdict(dict)
    all_pred_genera = defaultdict(dict)
    all_pred_full = {}
    sample_metrics = []
    
    print(f"Processing {args.samples} samples from {base_dir.absolute()}")
    print(f"Filters: min_kmers >= {args.min_kmers}, min_cov >= {args.min_cov}, min_reads >= {args.min_reads}")
    print()
    
    # Parse all files
    for i in range(1, args.samples + 1):
        sample_name = f"meta_sample_{i}"
        truth_file, report_file = find_files(base_dir, i)
        
        if not truth_file:
            print(f"  Warning: No truth file for sample {i}, skipping")
            continue
        if not report_file:
            print(f"  Warning: No report file for sample {i}, skipping")
            continue
        
        truth_genera = parse_taxonomic_profile(truth_file)
        pred_genera = parse_kraken_report(report_file)
        
        all_pred_full[sample_name] = pred_genera
        
        # Apply filters
        pred_filtered = apply_filter(
            pred_genera,
            min_reads=args.min_reads,
            min_kmers=args.min_kmers,
            min_cov=args.min_cov,
            min_taxReads=args.min_taxReads,
            max_dup=args.max_dup
        )
        
        for genus, abundance in truth_genera.items():
            all_truth_genera[genus][sample_name] = abundance
        for genus, data in pred_filtered.items():
            all_pred_genera[genus][sample_name] = data['reads']
        
        truth_set = set(truth_genera.keys())
        pred_set = set(pred_filtered.keys())
        
        binary_metrics = calculate_metrics_binary(truth_set, pred_set)
        abundance_metrics = calculate_metrics_abundance(truth_genera, {k: v['reads'] for k, v in pred_filtered.items()})
        
        sample_metrics.append({
            'Sample': sample_name,
            'Truth_Genera': len(truth_set),
            'Predicted_Genera': len(pred_set),
            **binary_metrics,
            **abundance_metrics
        })
        
        if i % 10 == 0:
            print(f"  Processed {i}/{args.samples} samples")
    
    print(f"\nProcessed {len(sample_metrics)} samples successfully")
    
    # Build matrices
    all_samples = [f"meta_sample_{i}" for i in range(1, args.samples + 1)]
    all_genera_union = sorted(set(all_truth_genera.keys()) | set(all_pred_genera.keys()))
    
    truth_matrix = pd.DataFrame(index=all_samples, columns=all_genera_union, dtype=float).fillna(0)
    for genus, samples_dict in all_truth_genera.items():
        for sample, abundance in samples_dict.items():
            if sample in truth_matrix.index:
                truth_matrix.loc[sample, genus] = abundance
    
    pred_matrix = pd.DataFrame(index=all_samples, columns=all_genera_union, dtype=float).fillna(0)
    for genus, samples_dict in all_pred_genera.items():
        for sample, reads in samples_dict.items():
            if sample in pred_matrix.index:
                pred_matrix.loc[sample, genus] = reads
    
    truth_pa = (truth_matrix > 0).astype(int)
    pred_pa = (pred_matrix > 0).astype(int)
    metrics_df = pd.DataFrame(sample_metrics)
    
    # Summary
    print("\n" + "="*70)
    print("SUMMARY STATISTICS")
    print("="*70)
    print(f"\nFilters applied: kmers >= {args.min_kmers}, cov >= {args.min_cov}")
    print(f"\nTotal unique genera in ground truth: {len(all_truth_genera)}")
    print(f"Total unique genera in predictions:  {len(all_pred_genera)}")
    print(f"Genera in both:                      {len(set(all_truth_genera.keys()) & set(all_pred_genera.keys()))}")
    print(f"\nBinary Metrics (mean +/- std across {len(sample_metrics)} samples):")
    print(f"  Precision (PPV):  {metrics_df['Precision_PPV'].mean():.4f} +/- {metrics_df['Precision_PPV'].std():.4f}")
    print(f"  Recall:           {metrics_df['Recall_Sensitivity'].mean():.4f} +/- {metrics_df['Recall_Sensitivity'].std():.4f}")
    print(f"  F1 Score:         {metrics_df['F1'].mean():.4f} +/- {metrics_df['F1'].std():.4f}")
    print(f"\nAbundance Metrics (mean +/- std):")
    print(f"  Bray-Curtis:      {metrics_df['Bray_Curtis_Dissimilarity'].mean():.4f} +/- {metrics_df['Bray_Curtis_Dissimilarity'].std():.4f}")
    print(f"  Correlation:      {metrics_df['Abundance_Correlation'].mean():.4f} +/- {metrics_df['Abundance_Correlation'].std():.4f}")
    
    # Save outputs
    output_prefix = args.output
    truth_matrix.to_csv(f"{output_prefix}_ground_truth_abundance.csv")
    pred_matrix.to_csv(f"{output_prefix}_krakenuniq_reads.csv")
    truth_pa.to_csv(f"{output_prefix}_ground_truth_presence_absence.csv")
    pred_pa.to_csv(f"{output_prefix}_krakenuniq_presence_absence.csv")
    metrics_df.to_csv(f"{output_prefix}_per_sample_metrics.csv", index=False)
    
    # Threshold analysis
    print("\n" + "="*70)
    print("THRESHOLD COMPARISON")
    print("="*70)
    
    threshold_results = []
    
    # Test key combinations
    test_configs = [
        {'name': 'No filter', 'min_kmers': 0, 'min_cov': 0},
        {'name': 'cov >= 0.001 only', 'min_kmers': 0, 'min_cov': 0.001},
        {'name': 'kmers >= 100 only', 'min_kmers': 100, 'min_cov': 0},
        {'name': 'OPTIMAL: kmers >= 100, cov >= 0.001', 'min_kmers': 100, 'min_cov': 0.001},
        {'name': 'kmers >= 500, cov >= 0.001', 'min_kmers': 500, 'min_cov': 0.001},
        {'name': 'kmers >= 100, cov >= 0.002', 'min_kmers': 100, 'min_cov': 0.002},
        {'name': 'High precision: cov >= 0.005', 'min_kmers': 0, 'min_cov': 0.005},
        {'name': 'reads >= 100, cov >= 0.001', 'min_reads': 100, 'min_kmers': 0, 'min_cov': 0.001},
    ]
    
    print(f"\n{'Configuration':<45} {'Precision':>10} {'Recall':>10} {'F1':>10}")
    print("-" * 80)
    
    for config in test_configs:
        name = config.pop('name')
        result = test_threshold(all_pred_full, truth_matrix, args.samples, **config)
        if result:
            threshold_results.append({'Config': name, **config, **result})
            print(f"{name:<45} {result['Precision']:>10.4f} {result['Recall']:>10.4f} {result['F1']:>10.4f}")
    
    threshold_df = pd.DataFrame(threshold_results)
    threshold_df.to_csv(f"{output_prefix}_threshold_comparison.csv", index=False)
    
    # Save Excel
    with pd.ExcelWriter(f"{output_prefix}_full_results.xlsx", engine='openpyxl') as writer:
        metrics_df.to_excel(writer, sheet_name='Per_Sample_Metrics', index=False)
        threshold_df.to_excel(writer, sheet_name='Threshold_Comparison', index=False)
        truth_matrix.to_excel(writer, sheet_name='Ground_Truth_Abundance')
        pred_matrix.to_excel(writer, sheet_name='KrakenUniq_Reads')
        truth_pa.to_excel(writer, sheet_name='Ground_Truth_PA')
        pred_pa.to_excel(writer, sheet_name='KrakenUniq_PA')
    
    print(f"\nOutputs saved:")
    print(f"  {output_prefix}_ground_truth_abundance.csv")
    print(f"  {output_prefix}_krakenuniq_reads.csv")
    print(f"  {output_prefix}_ground_truth_presence_absence.csv")
    print(f"  {output_prefix}_krakenuniq_presence_absence.csv")
    print(f"  {output_prefix}_per_sample_metrics.csv")
    print(f"  {output_prefix}_threshold_comparison.csv")
    print(f"  {output_prefix}_full_results.xlsx")


if __name__ == '__main__':
    main()
