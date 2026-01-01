#!/usr/bin/env python3
"""
Calculate Pearson and Spearman correlation between biological replicates
based on read counts in merged peaks.

This script:
1. Loads pre-computed count matrix (from bedtools coverage)
2. Normalizes to CPM (Counts Per Million) and log2 transforms
3. Calculates Pearson and Spearman correlation
4. Generates scatter plot visualization
"""

import argparse
import sys
import pandas as pd
import numpy as np
import scipy.stats as stats
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import matplotlib.pyplot as plt

# Check if running under Snakemake
try:
    snakemake
except NameError:
    snakemake = None


def calculate_correlation(counts_file, out_tsv, out_png):
    """
    Calculate correlation metrics from count matrix.
    
    Args:
        counts_file: TSV file with columns: chrom, start, end, name, rep1_count, rep2_count
        out_tsv: Output TSV file for correlation statistics
        out_png: Output PNG file for scatter plot
    """
    # Load count matrix
    try:
        df = pd.read_csv(
            counts_file, 
            sep='\t', 
            header=None, 
            names=['chrom', 'start', 'end', 'name', 'rep1', 'rep2']
        )
    except Exception as e:
        print(f"Error reading count file {counts_file}: {e}", file=sys.stderr)
        sys.exit(1)
    
    # Check for empty dataframe
    if df.empty:
        print("Warning: Empty count matrix", file=sys.stderr)
        # Create empty output files
        stats_df = pd.DataFrame({
            'metric': ['pearson_r', 'pearson_p', 'spearman_r', 'spearman_p'],
            'value': [np.nan, np.nan, np.nan, np.nan]
        })
        stats_df.to_csv(out_tsv, sep='\t', index=False)
        # Create empty plot
        fig, ax = plt.subplots(figsize=(6, 6))
        ax.text(0.5, 0.5, 'No data available', ha='center', va='center')
        plt.savefig(out_png, dpi=300, bbox_inches='tight')
        plt.close()
        return
    
    # Check for zero counts
    if df['rep1'].sum() == 0 or df['rep2'].sum() == 0:
        print("Warning: One or both replicates have zero total counts", file=sys.stderr)
        stats_df = pd.DataFrame({
            'metric': ['pearson_r', 'pearson_p', 'spearman_r', 'spearman_p'],
            'value': [np.nan, np.nan, np.nan, np.nan]
        })
        stats_df.to_csv(out_tsv, sep='\t', index=False)
        # Create empty plot
        fig, ax = plt.subplots(figsize=(6, 6))
        ax.text(0.5, 0.5, 'Zero counts detected', ha='center', va='center')
        plt.savefig(out_png, dpi=300, bbox_inches='tight')
        plt.close()
        return
    
    # 1. Normalization: Convert to CPM (Counts Per Million)
    total_reads_rep1 = df['rep1'].sum()
    total_reads_rep2 = df['rep2'].sum()
    
    cpm_rep1 = (df['rep1'] / total_reads_rep1) * 1e6
    cpm_rep2 = (df['rep2'] / total_reads_rep2) * 1e6
    
    # 2. Log transform: log2(CPM + 1)
    x = np.log2(cpm_rep1 + 1)
    y = np.log2(cpm_rep2 + 1)
    
    # 3. Calculate correlation statistics
    try:
        pearson_r, pearson_p = stats.pearsonr(x, y)
    except Exception as e:
        print(f"Error calculating Pearson correlation: {e}", file=sys.stderr)
        pearson_r, pearson_p = np.nan, np.nan
    
    try:
        spearman_r, spearman_p = stats.spearmanr(x, y)
    except Exception as e:
        print(f"Error calculating Spearman correlation: {e}", file=sys.stderr)
        spearman_r, spearman_p = np.nan, np.nan
    
    # 4. Save statistics to TSV
    stats_df = pd.DataFrame({
        'metric': ['pearson_r', 'pearson_p', 'spearman_r', 'spearman_p'],
        'value': [pearson_r, pearson_p, spearman_r, spearman_p]
    })
    stats_df.to_csv(out_tsv, sep='\t', index=False)
    
    # 5. Generate scatter plot
    fig, ax = plt.subplots(figsize=(6, 6))
    
    # Use hexbin for density visualization
    hb = ax.hexbin(x, y, gridsize=50, cmap='Blues', mincnt=1)
    plt.colorbar(hb, ax=ax, label='Density')
    
    # Add identity line (y = x)
    max_val = max(x.max(), y.max())
    ax.plot([0, max_val], [0, max_val], 'r--', linewidth=1.5, label='Identity line')
    
    # Labels and title
    ax.set_xlabel('Rep1 (log2 CPM + 1)', fontsize=12)
    ax.set_ylabel('Rep2 (log2 CPM + 1)', fontsize=12)
    ax.set_title(f'Replicate Correlation\nPearson R = {pearson_r:.4f}', fontsize=12)
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(out_png, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Correlation calculated: Pearson R = {pearson_r:.4f}, Spearman R = {spearman_r:.4f}", file=sys.stderr)


def main():
    # If running under Snakemake, use snakemake object
    if snakemake is not None:
        counts_file = snakemake.input.counts
        out_tsv = str(snakemake.output.tsv)
        out_png = str(snakemake.output.png)
    else:
        # Standalone mode: use argparse
        parser = argparse.ArgumentParser(
            description='Calculate correlation between biological replicates'
        )
        parser.add_argument(
            '--input', '-i',
            required=True,
            help='Input count matrix TSV (chrom, start, end, name, rep1_count, rep2_count)'
        )
        parser.add_argument(
            '--out-tsv', '-o',
            required=True,
            help='Output TSV file for correlation statistics'
        )
        parser.add_argument(
            '--out-png', '-p',
            required=True,
            help='Output PNG file for scatter plot'
        )
        
        args = parser.parse_args()
        counts_file = args.input
        out_tsv = args.out_tsv
        out_png = args.out_png
    
    calculate_correlation(counts_file, out_tsv, out_png)


if __name__ == '__main__':
    main()


