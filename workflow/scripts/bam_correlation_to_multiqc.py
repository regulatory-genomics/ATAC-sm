#!/usr/bin/env python3
"""
Convert bamReproducibility correlation matrices to MultiQC custom content format.
Calculates mean of off-diagonal values as the final metric for each replicate group.
"""

import argparse
import sys
from pathlib import Path
import yaml
import numpy as np


def parse_args():
    parser = argparse.ArgumentParser(
        description='Convert BAM correlation matrices to MultiQC format'
    )
    parser.add_argument(
        'matrix_files',
        nargs='+',
        help='Input correlation matrix files (*_global_rep_cor.txt)'
    )
    parser.add_argument(
        '--output-yaml',
        required=True,
        help='Output YAML file for MultiQC custom content'
    )
    parser.add_argument(
        '--output-stats',
        required=True,
        help='Output TSV file for MultiQC general statistics'
    )
    return parser.parse_args()


def parse_correlation_matrix(matrix_file):
    """
    Parse deepTools plotCorrelation output matrix.
    
    Format:
    	sample1	sample2	sample3
    sample1	1.0	0.95	0.93
    sample2	0.95	1.0	0.94
    sample3	0.93	0.94	1.0
    
    Or with quotes:
    	'test2.filtered.bam'	'test.filtered.bam'
    'test2.filtered.bam'	1.0000	1.0000
    'test.filtered.bam'	1.0000	1.0000
    """
    def strip_quotes(s):
        """Strip single or double quotes from string"""
        s = s.strip()
        if (s.startswith("'") and s.endswith("'")) or (s.startswith('"') and s.endswith('"')):
            return s[1:-1]
        return s
    
    def try_float(v):
        """Try to convert to float, return None if fails"""
        try:
            return float(v)
        except (ValueError, TypeError):
            return None
    
    with open(matrix_file, 'r') as f:
        lines = f.readlines()
    
    # Filter out comment lines (starting with #) and empty lines
    lines = [line for line in lines if line.strip() and not line.strip().startswith('#')]
    
    if len(lines) < 2:
        raise ValueError(f"Matrix file {matrix_file} has fewer than 2 non-comment lines")
    
    # Parse header (sample names) - first non-comment line
    header = lines[0].strip().split('\t')
    # Skip first empty column if present, strip quotes from sample names
    sample_names = [strip_quotes(h) for h in header[1:] if h.strip()]
    
    # Parse matrix values
    matrix = []
    row_names = []
    for line in lines[1:]:
        parts = line.strip().split('\t')
        if len(parts) < 2:
            continue
        
        # First column is row name (sample name), strip quotes
        row_name = strip_quotes(parts[0])
        if not row_name:
            continue
        row_names.append(row_name)
        
        # Remaining columns are correlation values (skip first column which is row name)
        values = []
        for v in parts[1:]:
            v = v.strip()
            if not v:
                continue
            # Try to convert to float
            val = try_float(v)
            if val is not None:
                values.append(val)
            else:
                # If it's a quoted string, try stripping quotes and parsing
                stripped = strip_quotes(v)
                val = try_float(stripped)
                if val is not None:
                    values.append(val)
                # Otherwise skip (might be a sample name in wrong position)
        
        if len(values) > 0:
            matrix.append(values)
    
    if len(matrix) == 0:
        raise ValueError(f"Could not parse any numeric values from {matrix_file}")
    
    matrix = np.array(matrix)
    
    # Verify dimensions match
    if matrix.shape[0] != len(sample_names):
        # Adjust sample_names to match matrix dimensions
        if len(sample_names) > matrix.shape[0]:
            sample_names = sample_names[:matrix.shape[0]]
        elif len(sample_names) < matrix.shape[0]:
            # Use row_names if available
            if len(row_names) == matrix.shape[0]:
                sample_names = row_names
            else:
                sample_names = [f"sample_{i+1}" for i in range(matrix.shape[0])]
    
    return sample_names, matrix


def calculate_off_diagonal_mean(matrix):
    """
    Calculate mean of off-diagonal elements.
    This represents the average correlation between different replicates.
    """
    n = matrix.shape[0]
    if n < 2:
        return None
    
    # Create mask for off-diagonal elements
    mask = ~np.eye(n, dtype=bool)
    
    # Get off-diagonal values
    off_diagonal = matrix[mask]
    
    # Calculate mean
    mean_corr = np.mean(off_diagonal)
    
    return mean_corr


def calculate_row_off_diagonal_means(matrix, sample_names):
    """
    Calculate mean of off-diagonal values for each row (sample).
    This shows how well each sample correlates with others.
    """
    n = matrix.shape[0]
    row_means = {}
    
    for i in range(n):
        # Get all values in row i except diagonal
        row = matrix[i, :]
        off_diag_values = np.concatenate([row[:i], row[i+1:]])
        
        if len(off_diag_values) > 0:
            row_means[sample_names[i]] = np.mean(off_diag_values)
        else:
            row_means[sample_names[i]] = None
    
    return row_means


def get_correlation_quality(mean_corr):
    """Determine quality level based on correlation value"""
    if mean_corr is None:
        return 'unknown'
    elif mean_corr >= 0.95:
        return 'excellent'
    elif mean_corr >= 0.90:
        return 'good'
    elif mean_corr >= 0.80:
        return 'acceptable'
    else:
        return 'poor'


def main():
    args = parse_args()
    
    all_data = {}
    general_stats = {}
    detailed_data = {}
    
    for matrix_file in args.matrix_files:
        matrix_path = Path(matrix_file)
        if not matrix_path.exists():
            print(f"Warning: {matrix_file} does not exist, skipping", file=sys.stderr)
            continue
        
        # Extract group name from filename
        # Assuming format: {group}_global_rep_cor.txt
        group_name = matrix_path.stem.replace('_global_rep_cor', '')
        
        try:
            # Parse correlation matrix
            sample_names, matrix = parse_correlation_matrix(matrix_file)
            
            # Verify matrix dimensions
            if matrix.shape[0] != matrix.shape[1]:
                print(f"Warning: Matrix is not square ({matrix.shape[0]}x{matrix.shape[1]}) for {group_name}", file=sys.stderr)
            
            if len(sample_names) != matrix.shape[0]:
                print(f"Warning: Number of sample names ({len(sample_names)}) doesn't match matrix rows ({matrix.shape[0]}) for {group_name}", file=sys.stderr)
            
            # Calculate overall mean of off-diagonal values
            mean_corr = calculate_off_diagonal_mean(matrix)
            
            # Calculate per-sample means
            row_means = calculate_row_off_diagonal_means(matrix, sample_names)
            
            # Determine quality
            quality = get_correlation_quality(mean_corr)
            
            # Store for general statistics
            if mean_corr is not None:
                general_stats[group_name] = {
                    'bam_correlation': mean_corr,
                    'num_replicates': len(sample_names),
                    'quality': quality
                }
                
                # Store detailed data
                all_data[group_name] = {
                    'Mean_Correlation': f"{mean_corr:.4f}",
                    'Quality': quality,
                    'Num_Replicates': len(sample_names),
                    'Min_Correlation': f"{np.min(matrix[~np.eye(len(sample_names), dtype=bool)]):.4f}",
                    'Max_Correlation': f"{np.max(matrix[~np.eye(len(sample_names), dtype=bool)]):.4f}",
                }
                
                # Add per-sample correlations
                for sample, mean_val in row_means.items():
                    if mean_val is not None:
                        detailed_data[f"{group_name}_{sample}"] = {
                            'Group': group_name,
                            'Sample': sample,
                            'Mean_Correlation': f"{mean_val:.4f}",
                        }
            
            print(f"Processed {group_name}: mean_corr={mean_corr:.4f}, "
                  f"n_replicates={len(sample_names)}, quality={quality}", 
                  file=sys.stderr)
            
        except Exception as e:
            print(f"Error processing {matrix_file}: {e}", file=sys.stderr)
            continue
    
    if not all_data:
        print("Error: No valid matrix files found", file=sys.stderr)
        sys.exit(1)
    
    # Create MultiQC custom content YAML for group-level summary
    multiqc_data_groups = {
        'id': 'bam_correlation_groups',
        'section_name': 'BAM Correlation - Replicate Groups',
        'description': 'Spearman correlation between BAM files within replicate groups. Mean correlation is calculated from off-diagonal values (excluding self-correlation).',
        'plot_type': 'table',
        'pconfig': {
            'id': 'bam_correlation_table',
            'title': 'BAM Correlation Metrics',
            'save_file': True,
            'col1_header': 'Replicate Group',
        },
        'headers': {
            'Mean_Correlation': {
                'title': 'Mean Correlation',
                'description': 'Mean of off-diagonal Spearman correlation values',
                'min': 0,
                'max': 1,
                'format': '{:.4f}',
                'scale': 'RdYlGn',
            },
            'Quality': {
                'title': 'Quality',
                'description': 'Quality assessment: excellent (≥0.95), good (≥0.90), acceptable (≥0.80), poor (<0.80)',
                'scale': False,
            },
            'Num_Replicates': {
                'title': 'N Replicates',
                'description': 'Number of replicates in the group',
                'format': '{:d}',
            },
            'Min_Correlation': {
                'title': 'Min Correlation',
                'description': 'Minimum pairwise correlation',
                'format': '{:s}',
                'hidden': True,
            },
            'Max_Correlation': {
                'title': 'Max Correlation',
                'description': 'Maximum pairwise correlation',
                'format': '{:s}',
                'hidden': True,
            },
        },
        'data': all_data,
    }
    
    # Create MultiQC custom content YAML for per-sample details
    multiqc_data_samples = {
        'id': 'bam_correlation_samples',
        'section_name': 'BAM Correlation - Individual Samples',
        'description': 'Mean correlation of each sample with other samples in its replicate group.',
        'plot_type': 'table',
        'pconfig': {
            'id': 'bam_correlation_sample_table',
            'title': 'Per-Sample Correlation',
            'save_file': True,
            'col1_header': 'Sample',
        },
        'headers': {
            'Group': {
                'title': 'Group',
                'description': 'Replicate group',
                'scale': False,
            },
            'Sample': {
                'title': 'Sample',
                'description': 'Sample name',
                'scale': False,
            },
            'Mean_Correlation': {
                'title': 'Mean Correlation',
                'description': 'Mean correlation with other samples in group',
                'min': 0,
                'max': 1,
                'format': '{:s}',
                'scale': 'RdYlGn',
            },
        },
        'data': detailed_data,
    }
    
    # Combine both sections
    multiqc_combined = [multiqc_data_groups, multiqc_data_samples]
    
    # Write MultiQC YAML
    output_yaml_path = Path(args.output_yaml)
    output_yaml_path.parent.mkdir(parents=True, exist_ok=True)
    
    with open(args.output_yaml, 'w') as f:
        yaml.dump(multiqc_combined, f, default_flow_style=False, sort_keys=False)
    
    print(f"Created MultiQC custom content: {args.output_yaml}", file=sys.stderr)
    
    # Create general statistics TSV (for main MultiQC table)
    output_stats_path = Path(args.output_stats)
    output_stats_path.parent.mkdir(parents=True, exist_ok=True)
    
    with open(args.output_stats, 'w') as f:
        # Header
        f.write('Sample\tBAM_Correlation\tCorrelation_Quality\tN_Replicates\n')
        
        # Data rows
        for group_name, stats in general_stats.items():
            f.write(f"{group_name}\t"
                   f"{stats['bam_correlation']:.4f}\t"
                   f"{stats['quality']}\t"
                   f"{stats['num_replicates']}\n")
    
    print(f"Created general statistics TSV: {args.output_stats}", file=sys.stderr)
    print(f"Processed {len(all_data)} replicate groups", file=sys.stderr)


if __name__ == '__main__':
    main()

