#!/usr/bin/env python3
"""
Reproducibility QC for ATAC-seq peaks following ENCODE standards.

Evaluates peak consistency and selects optimal peak set based on:
1. Self-consistency (pr1 vs pr2 for each replicate)
2. True replicate consistency (between biological replicates)
3. Pooled pseudo-replicate consistency (pooled-pr1 vs pooled-pr2)

ENCODE optimal peak selection:
- If Nt >= Np: use Nt (max of true replicate pairs)
- If Nt < Np: use Np (pooled pseudo-replicates)
"""

import argparse
import gzip
import json
import sys
from pathlib import Path

# Check if running under Snakemake
try:
    snakemake
except NameError:
    snakemake = None


def count_peaks(peak_file):
    """Count number of peaks in a narrowPeak file"""
    if not Path(peak_file).exists():
        return 0
    
    count = 0
    opener = gzip.open if peak_file.endswith('.gz') else open
    with opener(peak_file, 'rt') as f:
        for line in f:
            if line.strip() and not line.startswith('#'):
                count += 1
    return count


def calculate_rescue_ratio(n1, n2, np):
    """
    Calculate rescue ratio: max(N1, N2) / Np
    A rescue ratio > 2 indicates poor replicate consistency
    """
    if np == 0:
        return float('inf')
    return max(n1, n2) / np


def calculate_self_consistency_ratio(n1, n2):
    """
    Calculate self-consistency ratio: max(N1, N2) / min(N1, N2)
    Closer to 1 = better consistency
    """
    if min(n1, n2) == 0:
        return float('inf')
    return max(n1, n2) / min(n1, n2)


def main():
    # If running under Snakemake, use snakemake object
    if snakemake is not None:
        # Extract inputs from snakemake object
        # peaks_pr is a list, peak_ppr is a single file, true_rep_pairs is a list (not used in script)
        peaks_pr = list(snakemake.input.peaks_pr)
        peak_ppr = snakemake.input.peak_ppr
        prefix = snakemake.params.prefix
        output_qc = str(snakemake.output.qc_json)
        output_optimal = str(snakemake.output.optimal_peak)
    else:
        # Standalone mode: use argparse
        parser = argparse.ArgumentParser(description='Reproducibility QC for ATAC-seq peaks')
        parser.add_argument('--peaks-pr', nargs='+', required=True,
                            help='Self-consistency peaks (pr1 vs pr2 for each replicate)')
        parser.add_argument('--peak-ppr', required=False,
                            help='Pooled pseudo-replicate peak file (pooled-pr1 vs pooled-pr2)')
        parser.add_argument('--prefix', required=True,
                            help='Output prefix (e.g., "idr" or "overlap")')
        parser.add_argument('--output-qc', required=True,
                            help='Output QC JSON file')
        parser.add_argument('--output-optimal', required=True,
                            help='Output optimal peak file (symlink to best peak set)')
        
        args = parser.parse_args()
        peaks_pr = args.peaks_pr
        peak_ppr = args.peak_ppr
        prefix = args.prefix
        output_qc = args.output_qc
        output_optimal = args.output_optimal
    
    # Count peaks in self-consistency files (pr1 vs pr2 for each replicate)
    self_consistency_counts = []
    for peak_file in peaks_pr:
        count = count_peaks(peak_file)
        self_consistency_counts.append({
            'file': Path(peak_file).name,
            'num_peaks': count
        })
        print(f"Self-consistency: {Path(peak_file).name} = {count} peaks", file=sys.stderr)
    
    # Count peaks in pooled pseudo-replicates (if provided)
    if peak_ppr and Path(peak_ppr).exists():
        np_count = count_peaks(peak_ppr)
        print(f"Pooled PR: {Path(peak_ppr).name} = {np_count} peaks", file=sys.stderr)
    else:
        np_count = 0
        print("Pooled PR: Not available", file=sys.stderr)
    
    # Find maximum self-consistency peak count (this represents Nt in ENCODE terminology)
    if self_consistency_counts:
        max_self_consistency = max(self_consistency_counts, key=lambda x: x['num_peaks'])
        nt_count = max_self_consistency['num_peaks']
        optimal_file = None
        
        # Find the corresponding peak file
        for peak_file in peaks_pr:
            if Path(peak_file).name == max_self_consistency['file']:
                optimal_file = peak_file
                break
    else:
        nt_count = 0
        max_self_consistency = None
        optimal_file = None
    
    # Determine optimal peak set based on ENCODE criteria
    if nt_count >= np_count:
        # Use max of self-consistency (true replicate comparison)
        optimal_peak_file = optimal_file
        optimal_peak_type = 'self_consistency'
        optimal_peak_count = nt_count
    else:
        # Use pooled pseudo-replicates
        optimal_peak_file = peak_ppr if peak_ppr else optimal_file
        optimal_peak_type = 'pooled_pseudoreplicates'
        optimal_peak_count = np_count
    
    print(f"\nOptimal peak: {optimal_peak_type} with {optimal_peak_count} peaks", file=sys.stderr)
    
    # Calculate rescue ratio if we have multiple replicates
    rescue_ratio = None
    if len(self_consistency_counts) >= 2 and np_count > 0:
        # Sort to get top 2
        sorted_counts = sorted([x['num_peaks'] for x in self_consistency_counts], reverse=True)
        rescue_ratio = calculate_rescue_ratio(sorted_counts[0], sorted_counts[1], np_count)
    
    # Calculate self-consistency ratios
    self_consistency_ratios = []
    if len(self_consistency_counts) >= 2:
        counts = [x['num_peaks'] for x in self_consistency_counts]
        for i in range(len(counts)):
            for j in range(i+1, len(counts)):
                ratio = calculate_self_consistency_ratio(counts[i], counts[j])
                self_consistency_ratios.append({
                    'rep_i': i+1,
                    'rep_j': j+1,
                    'ratio': ratio
                })
    
    # Create QC output
    qc_output = {
        'prefix': prefix,
        'num_replicates': len(self_consistency_counts),
        'self_consistency': self_consistency_counts,
        'pooled_pseudoreplicates': {
            'file': Path(peak_ppr).name if peak_ppr else None,
            'num_peaks': np_count
        },
        'optimal_peak': {
            'type': optimal_peak_type,
            'file': Path(optimal_peak_file).name if optimal_peak_file else None,
            'num_peaks': optimal_peak_count
        },
        'Nt': nt_count,
        'Np': np_count,
        'rescue_ratio': rescue_ratio,
        'self_consistency_ratios': self_consistency_ratios,
        'reproducibility': 'PASS' if (rescue_ratio is None or rescue_ratio <= 2) else 'FAIL'
    }
    
    # Write QC JSON
    with open(output_qc, 'w') as f:
        json.dump(qc_output, f, indent=2)
    
    # Create symlink to optimal peak file
    if optimal_peak_file and Path(optimal_peak_file).exists():
        output_path = Path(output_optimal)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        # Remove existing symlink/file
        if output_path.exists() or output_path.is_symlink():
            output_path.unlink()
        
        # Create symlink
        output_path.symlink_to(Path(optimal_peak_file).absolute())
        print(f"Created symlink: {output_optimal} -> {optimal_peak_file}", file=sys.stderr)
    else:
        # Create empty file if no optimal peak available
        Path(output_optimal).touch()
        print(f"Warning: No optimal peak file available, created empty file", file=sys.stderr)
    
    print(f"\nReproducibility: {qc_output['reproducibility']}", file=sys.stderr)
    if rescue_ratio is not None:
        print(f"Rescue ratio: {rescue_ratio:.2f} ({'PASS' if rescue_ratio <= 2 else 'FAIL'})", file=sys.stderr)


if __name__ == '__main__':
    main()

