#!/usr/bin/env python3
"""
Aggregate correlation results from all replicate pairs in a group.

This script collects correlation statistics (Pearson and Spearman) from
individual pair-wise correlation TSV files and aggregates them into a
single summary file.
"""

import pandas as pd
import os
import sys

# Check if running under Snakemake
try:
    snakemake
except NameError:
    snakemake = None


def aggregate_correlations(pairs, out_file):
    """
    Aggregate correlation statistics from multiple pair files.
    
    Args:
        pairs: List of tuples (rep1, rep2, file_path)
        out_file: Output TSV file path
    """
    results = []
    
    for rep1, rep2, path in pairs:
        if os.path.exists(path):
            try:
                # Read correlation statistics
                df = pd.read_csv(path, sep='\t')
                # Convert to dictionary for easier access
                stats_dict = dict(zip(df['metric'], df['value']))
                
                results.append({
                    'rep1': rep1,
                    'rep2': rep2,
                    'pearson_r': stats_dict.get('pearson_r', None),
                    'pearson_p': stats_dict.get('pearson_p', None),
                    'spearman_r': stats_dict.get('spearman_r', None),
                    'spearman_p': stats_dict.get('spearman_p', None),
                })
            except Exception as e:
                print(f"Error reading {path}: {e}", file=sys.stderr)
        else:
            print(f"Warning: Correlation file not found: {path}", file=sys.stderr)
    
    # Create summary dataframe
    if results:
        summary_df = pd.DataFrame(results)
    else:
        # Create empty dataframe with correct columns
        summary_df = pd.DataFrame(columns=['rep1', 'rep2', 'pearson_r', 'pearson_p', 'spearman_r', 'spearman_p'])
    
    # Save summary
    os.makedirs(os.path.dirname(out_file), exist_ok=True)
    summary_df.to_csv(out_file, sep='\t', index=False)
    
    print(f"Aggregated {len(results)} correlation results to {out_file}", file=sys.stderr)


def main():
    # If running under Snakemake, use snakemake object
    if snakemake is not None:
        pairs = snakemake.params.pairs
        out_file = str(snakemake.output.summary_tsv)
    else:
        # Standalone mode: use argparse
        import argparse
        parser = argparse.ArgumentParser(
            description='Aggregate correlation results from replicate pairs'
        )
        parser.add_argument(
            '--pairs',
            nargs='+',
            required=True,
            help='List of tuples: rep1,rep2,file_path (space-separated)'
        )
        parser.add_argument(
            '--output', '-o',
            required=True,
            help='Output TSV file for aggregated correlations'
        )
        
        args = parser.parse_args()
        # Parse pairs from command line (format: rep1,rep2,path rep1,rep2,path ...)
        pairs = []
        for pair_str in args.pairs:
            parts = pair_str.split(',')
            if len(parts) == 3:
                pairs.append((parts[0], parts[1], parts[2]))
            else:
                print(f"Warning: Invalid pair format: {pair_str}", file=sys.stderr)
        out_file = args.output
    
    aggregate_correlations(pairs, out_file)


if __name__ == '__main__':
    main()

