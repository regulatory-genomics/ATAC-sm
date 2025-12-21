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


def calculate_rescue_ratio(np, nt):
    """
    Calculate rescue ratio: max(Np, Nt) / min(Np, Nt)
    A rescue ratio > 2 indicates poor consistency
    """
    if min(np, nt) == 0:
        return float('inf')
    return float(max(np, nt)) / float(min(np, nt))


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
        peaks_pr = list(snakemake.input.peaks_pr)  # Self-consistency: N1, N2, ...
        peak_ppr = snakemake.input.peak_ppr  # Pooled pseudo-replicates: Np
        peaks_true_reps = list(snakemake.input.true_rep_pairs) if hasattr(snakemake.input, 'true_rep_pairs') else []  # True replicates: Nt
        prefix = snakemake.params.prefix
        output_qc = str(snakemake.output.qc_json)
        output_optimal = str(snakemake.output.optimal_peak)
    else:
        # Standalone mode: use argparse
        parser = argparse.ArgumentParser(description='Reproducibility QC for ATAC-seq peaks')
        parser.add_argument('--peaks-pr', nargs='+', required=True,
                            help='Self-consistency peaks (pr1 vs pr2 for each replicate) - N1, N2, ...')
        parser.add_argument('--peaks-true-reps', nargs='*', default=[],
                            help='True replicate peaks (rep1_vs_rep2, etc.) - for Nt calculation')
        parser.add_argument('--peak-ppr', required=False,
                            help='Pooled pseudo-replicate peak file (pooled-pr1 vs pooled-pr2) - Np')
        parser.add_argument('--prefix', required=True,
                            help='Output prefix')
        parser.add_argument('--output-qc', required=True,
                            help='Output QC JSON file')
        parser.add_argument('--output-optimal', required=True,
                            help='Output optimal peak file (symlink to best peak set)')
        
        args = parser.parse_args()
        peaks_pr = args.peaks_pr
        peaks_true_reps = args.peaks_true_reps
        peak_ppr = args.peak_ppr
        prefix = args.prefix
        output_qc = args.output_qc
        output_optimal = args.output_optimal
    
    # ========================================================================
    # ENCODE Definitions:
    # N1, N2, ... = Self-consistency for each replicate (rep1-pr1 vs rep1-pr2)
    # Nt = True replicate consistency (max of rep1_vs_rep2, rep1_vs_rep3, etc.)
    # Np = Pooled pseudo-replicate consistency (pooled-pr1 vs pooled-pr2)
    # ========================================================================
    
    # 1. Count N1, N2, ... (Self-consistency for each replicate)
    self_consistency_counts = []
    n_self = []
    for peak_file in peaks_pr:
        count = count_peaks(peak_file)
        self_consistency_counts.append({
            'file': Path(peak_file).name,
            'num_peaks': count
        })
        n_self.append(count)
        print(f"Self-consistency: {Path(peak_file).name} = {count} peaks", file=sys.stderr)
    
    # 2. Count Nt (True replicate consistency - max of all pairwise comparisons)
    true_replicate_counts = []
    n_true = []
    for peak_file in peaks_true_reps:
        count = count_peaks(peak_file)
        true_replicate_counts.append({
            'file': Path(peak_file).name,
            'num_peaks': count
        })
        n_true.append(count)
        print(f"True replicate: {Path(peak_file).name} = {count} peaks", file=sys.stderr)
    
    nt_count = max(n_true) if n_true else 0
    
    # Find the conservative peak file (the one with Nt peaks)
    conservative_file = None
    if nt_count > 0:
        nt_idx = n_true.index(nt_count)
        conservative_file = peaks_true_reps[nt_idx]
        print(f"Conservative set (Nt={nt_count}): {Path(conservative_file).name}", file=sys.stderr)
    
    # 3. Count Np (Pooled pseudo-replicate consistency)
    if peak_ppr and Path(peak_ppr).exists():
        np_count = count_peaks(peak_ppr)
        print(f"Pooled pseudo-replicates (Np): {Path(peak_ppr).name} = {np_count} peaks", file=sys.stderr)
    else:
        np_count = 0
        print("Pooled PR (Np): Not available", file=sys.stderr)
    
    # 4. Determine optimal peak set based on ENCODE criteria
    if nt_count >= np_count:
        # Use conservative set (true replicates)
        optimal_peak_file = conservative_file
        optimal_peak_type = 'conservative_set'
        optimal_peak_count = nt_count
    else:
        # Use pooled pseudo-replicates
        optimal_peak_file = peak_ppr if peak_ppr else conservative_file
        optimal_peak_type = 'pooled_pseudoreplicates'
        optimal_peak_count = np_count
    
    print(f"\nOptimal peak: {optimal_peak_type} with {optimal_peak_count} peaks", file=sys.stderr)
    
    # 5. Calculate rescue ratio: max(Np, Nt) / min(Np, Nt)
    rescue_ratio = None
    if np_count > 0 and nt_count > 0:
        rescue_ratio = float(max(np_count, nt_count)) / float(min(np_count, nt_count))
    
    # 6. Calculate self-consistency ratio: max(N1, N2, ...) / min(N1, N2, ...)
    self_consistency_ratio = None
    if len(n_self) >= 2:
        self_consistency_ratio = float(max(n_self)) / float(min(n_self))
    
    # Determine reproducibility status
    reproducibility = 'pass'
    if rescue_ratio and self_consistency_ratio:
        if rescue_ratio > 2.0 or self_consistency_ratio > 2.0:
            reproducibility = 'borderline'
        if rescue_ratio > 2.0 and self_consistency_ratio > 2.0:
            reproducibility = 'fail'
    elif rescue_ratio and rescue_ratio > 2.0:
        reproducibility = 'borderline'
    elif self_consistency_ratio and self_consistency_ratio > 2.0:
        reproducibility = 'borderline'
    
    # Create QC output
    qc_output = {
        'prefix': prefix,
        'num_replicates': len(self_consistency_counts),
        'self_consistency': self_consistency_counts,
        'N_self': n_self,
        'true_replicates': true_replicate_counts,
        'N_true': n_true,
        'pooled_pseudoreplicates': {
            'file': Path(peak_ppr).name if peak_ppr else None,
            'num_peaks': np_count
        },
        'optimal_peak': {
            'type': optimal_peak_type,
            'file': Path(optimal_peak_file).name if optimal_peak_file else None,
            'num_peaks': optimal_peak_count
        },
        'conservative_peak': {
            'file': Path(conservative_file).name if conservative_file else None,
            'num_peaks': nt_count
        },
        'Nt': nt_count,
        'Np': np_count,
        'N_optimal': optimal_peak_count,
        'N_conservative': nt_count,
        'rescue_ratio': rescue_ratio,
        'self_consistency_ratio': self_consistency_ratio,
        'reproducibility': reproducibility
    }
    
    # Write QC JSON (ensure directory exists)
    output_qc_path = Path(output_qc)
    output_qc_path.parent.mkdir(parents=True, exist_ok=True)
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
    
    print(f"\nReproducibility: {reproducibility}", file=sys.stderr)
    if rescue_ratio is not None:
        print(f"Rescue ratio: {rescue_ratio:.2f} ({'pass' if rescue_ratio <= 2 else 'fail'})", file=sys.stderr)
    if self_consistency_ratio is not None:
        print(f"Self-consistency ratio: {self_consistency_ratio:.2f} ({'pass' if self_consistency_ratio <= 2 else 'fail'})", file=sys.stderr)


if __name__ == '__main__':
    main()

