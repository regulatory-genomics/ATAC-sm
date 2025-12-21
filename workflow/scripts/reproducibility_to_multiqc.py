#!/usr/bin/env python3
"""
Convert reproducibility QC JSON files to MultiQC custom content format.
Generates both a general statistics table and a detailed reproducibility metrics table.
"""

import json
import sys
import argparse
from pathlib import Path
import yaml


def parse_args():
    parser = argparse.ArgumentParser(
        description='Convert reproducibility QC JSON to MultiQC format'
    )
    parser.add_argument(
        'json_files',
        nargs='+',
        help='Input reproducibility QC JSON files'
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


def format_ratio(value):
    """Format ratio value for display"""
    if value is None:
        return 'N/A'
    if value == float('inf'):
        return 'inf'
    return f"{value:.3f}"


def get_reproducibility_color(status):
    """Get color based on reproducibility status"""
    colors = {
        'pass': '#5cb85c',      # green
        'borderline': '#f0ad4e',  # orange
        'fail': '#d9534f'         # red
    }
    return colors.get(status, '#999999')


def main():
    args = parse_args()
    
    # Collect data from all JSON files
    all_data = {}
    general_stats = {}
    
    for json_file in args.json_files:
        json_path = Path(json_file)
        if not json_path.exists():
            print(f"Warning: {json_file} does not exist, skipping", file=sys.stderr)
            continue
        
        with open(json_file, 'r') as f:
            data = json.load(f)
        
        # Extract sample/group name from filename
        # Assuming filename format: {group}.reproducibility.qc.json
        sample_name = json_path.stem.replace('.reproducibility.qc', '')
        
        # Extract metrics
        n_self_consistency = data.get('self_consistency', [])
        np_count = data.get('Np', 0)
        nt_count = data.get('Nt', 0)
        rescue_ratio = data.get('rescue_ratio')
        self_consistency_ratio = data.get('self_consistency_ratio')
        reproducibility = data.get('reproducibility', 'unknown')
        n_optimal = data.get('N_optimal', 0)
        
        # Store detailed data for custom content table
        all_data[sample_name] = {
            'Reproducibility': reproducibility,
            'N1': n_self_consistency[0]['num_peaks'] if len(n_self_consistency) > 0 else 0,
            'N2': n_self_consistency[1]['num_peaks'] if len(n_self_consistency) > 1 else 0,
            'Np': np_count,
            'Nt': nt_count,
            'Rescue_Ratio': format_ratio(rescue_ratio),
            'Self_Consistency_Ratio': format_ratio(self_consistency_ratio),
            'N_Optimal': n_optimal,
        }
        
        # Store data for general statistics (shown in MultiQC overview)
        general_stats[sample_name] = {
            'reproducibility': reproducibility,
            'n_optimal_peaks': n_optimal,
            'rescue_ratio': rescue_ratio if rescue_ratio and rescue_ratio != float('inf') else None,
            'self_consistency_ratio': self_consistency_ratio if self_consistency_ratio and self_consistency_ratio != float('inf') else None,
        }
    
    if not all_data:
        print("Error: No valid JSON files found", file=sys.stderr)
        sys.exit(1)
    
    # Create MultiQC custom content YAML
    multiqc_data = {
        'id': 'reproducibility_metrics',
        'section_name': 'ATAC-seq Reproducibility',
        'description': 'IDR-based reproducibility metrics for ATAC-seq peak calling. N1/N2: self-consistency peaks, Np: pooled pseudoreplicate peaks, Nt: true replicate consistency peaks.',
        'plot_type': 'table',
        'pconfig': {
            'id': 'reproducibility_table',
            'title': 'Reproducibility Metrics',
            'save_file': True,
            'col1_header': 'Sample Group',
        },
        'headers': {
            'Reproducibility': {
                'title': 'Status',
                'description': 'Overall reproducibility status',
                'scale': False,
                'format': '{:s}',
            },
            'N1': {
                'title': 'N1',
                'description': 'Self-consistency: Rep1-pr1 vs Rep1-pr2',
                'format': '{:,.0f}',
                'scale': 'Greens',
                'shared_key': 'peak_count',
            },
            'N2': {
                'title': 'N2',
                'description': 'Self-consistency: Rep2-pr1 vs Rep2-pr2',
                'format': '{:,.0f}',
                'scale': 'Greens',
                'shared_key': 'peak_count',
            },
            'Np': {
                'title': 'Np',
                'description': 'Pooled pseudoreplicates consistency',
                'format': '{:,.0f}',
                'scale': 'Blues',
                'shared_key': 'peak_count',
            },
            'Nt': {
                'title': 'Nt',
                'description': 'True replicate consistency',
                'format': '{:,.0f}',
                'scale': 'Purples',
                'shared_key': 'peak_count',
            },
            'Rescue_Ratio': {
                'title': 'Rescue Ratio',
                'description': 'max(Np,Nt) / min(Np,Nt). Should be ≤2.0',
                'format': '{:s}',
                'scale': False,
            },
            'Self_Consistency_Ratio': {
                'title': 'Self-Consist. Ratio',
                'description': 'max(N1,N2) / min(N1,N2). Should be ≤2.0',
                'format': '{:s}',
                'scale': False,
            },
            'N_Optimal': {
                'title': 'N Optimal',
                'description': 'Number of peaks in optimal set',
                'format': '{:,.0f}',
                'scale': 'RdYlGn',
                'shared_key': 'peak_count',
            },
        },
        'data': all_data,
    }
    
    # Write MultiQC YAML
    output_yaml_path = Path(args.output_yaml)
    output_yaml_path.parent.mkdir(parents=True, exist_ok=True)
    
    with open(args.output_yaml, 'w') as f:
        yaml.dump(multiqc_data, f, default_flow_style=False, sort_keys=False)
    
    print(f"Created MultiQC custom content: {args.output_yaml}", file=sys.stderr)
    
    # Create general statistics TSV (for main MultiQC table)
    output_stats_path = Path(args.output_stats)
    output_stats_path.parent.mkdir(parents=True, exist_ok=True)
    
    with open(args.output_stats, 'w') as f:
        # Header
        f.write('Sample\tReproducibility\tN_Optimal_Peaks\tRescue_Ratio\tSelf_Consistency_Ratio\n')
        
        # Data rows
        for sample_name, stats in general_stats.items():
            rescue = format_ratio(stats['rescue_ratio']) if stats['rescue_ratio'] is not None else 'N/A'
            self_consist = format_ratio(stats['self_consistency_ratio']) if stats['self_consistency_ratio'] is not None else 'N/A'
            
            f.write(f"{sample_name}\t"
                   f"{stats['reproducibility']}\t"
                   f"{stats['n_optimal_peaks']}\t"
                   f"{rescue}\t"
                   f"{self_consist}\n")
    
    print(f"Created general statistics TSV: {args.output_stats}", file=sys.stderr)
    print(f"Processed {len(all_data)} samples", file=sys.stderr)


if __name__ == '__main__':
    main()

