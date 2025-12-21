#!/usr/bin/env python

"""
MultiQC module to parse ATAC-seq pipeline stats
"""

from __future__ import print_function
from collections import OrderedDict
import logging
import pandas as pd
import os
import csv
import numpy as np

from multiqc import config
from multiqc.plots import linegraph, table, scatter
from multiqc.modules.base_module import BaseMultiqcModule

# Initialise the logger
log = logging.getLogger('multiqc')


class MultiqcModule(BaseMultiqcModule):
    """
    atacseq module class
    """

    def __init__(self):
        # Halt execution if we've disabled the plugin
        if config.kwargs.get('disable_atacseq_report', True):
            return None

        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='The ATAC-seq Pipeline', anchor='atacseq',
                                            href='https://github.com/regulatory-genomics/ATAC-sm',
                                            info="processes, quantifies and annotates ATAC-seq data.")
        log.info('Initialized atacseq module')
        
        # Parse ATAC-seq stats for each sample
        self.atacseq_data = dict()
        for f in self.find_log_files(sp_key='atacseq'):
            self.atacseq_data[f['s_name']] = self.parse_atacseq_stats(f['f'])
        log.info('Found stats file for {} ATAC-seq samples'.format(len(self.atacseq_data)))

        # Raise the not found warning
        if len(self.atacseq_data) == 0:
            raise UserWarning

        # Parse TSS for each sample
        self.atacseq_tss_data = dict()
        for f in self.find_log_files(sp_key='atacseq/tss'):
            my_sample_name = f['s_name'].replace('_TSS','')
            self.atacseq_tss_data[f['s_name']], self.atacseq_data[my_sample_name]['tss_max'] = self.parse_atacseq_tss(f['f'])
        log.info('Found TSS file for {} ATAC-seq samples'.format(len(self.atacseq_tss_data)))
        
        # Parse align_stats for each sample (contains mitochondrial_fraction)
        align_stats_files = list(self.find_log_files(sp_key='atacseq/align_stats'))
        for f in align_stats_files:
            # Extract sample name from filename (e.g., "sample.align.stats.tsv" -> "sample")
            s_name = f['s_name']
            if s_name.endswith('.align.stats'):
                sample_name = s_name.replace('.align.stats', '')
            elif s_name.endswith('.align'):
                sample_name = s_name.replace('.align', '')
            else:
                sample_name = s_name
            
            align_stats_data = self.parse_align_stats(f['f'])
            # Merge align_stats data into atacseq_data
            if sample_name not in self.atacseq_data:
                self.atacseq_data[sample_name] = {}
            self.atacseq_data[sample_name].update(align_stats_data)
        log.info('Found align_stats file for {} ATAC-seq samples'.format(len(align_stats_files)))
        
        # Parse Samblaster logs for n_tot, n_nondups, and NRF
        # Samblaster logs contain: "Marked X of Y (Z%) read ids as duplicates"
        # n_dups = X, n_tot = Y, n_nondups = Y - X
        samblaster_files = list(self.find_log_files('samblaster'))
        for f in samblaster_files:
            sample_name = f['s_name']
            # Parse the log file
            for line in f['f'].splitlines():
                if 'Marked' in line and 'read ids as duplicates' in line:
                    # Example: "samblaster: Marked 3894 of 7458 (52.21%) read ids as duplicates..."
                    import re
                    match = re.search(r'Marked\s+(\d+)\s+of\s+(\d+)', line)
                    if match:
                        n_dups = int(match.group(1))
                        n_tot = int(match.group(2))
                        n_nondups = n_tot - n_dups
                        
                        if sample_name not in self.atacseq_data:
                            self.atacseq_data[sample_name] = {}
                        
                        self.atacseq_data[sample_name]['n_tot'] = n_tot
                        self.atacseq_data[sample_name]['n_nondups'] = n_nondups
                        self.atacseq_data[sample_name]['n_dups'] = n_dups
                        
                        # Calculate NRF (Non-Redundant Fraction) = 1 - duplication_rate = n_nondups / n_tot
                        if n_tot > 0:
                            nrf = float(n_nondups) / float(n_tot)
                            self.atacseq_data[sample_name]['nrf'] = nrf
                        else:
                            self.atacseq_data[sample_name]['nrf'] = 0.0
                        break
        log.info('Found Samblaster log files for {} samples'.format(len(samblaster_files)))
        
        # Parse prealign stats files for percent_filtered
        # Prealign stats are per sample_run, need to map to sample_name
        prealign_stats_files = list(self.find_log_files(sp_key='atacseq/prealign_stats'))
        if prealign_stats_files:
            # Load annotation to map sample_run to sample_name
            try:
                import csv
                sample_annotation_path = config.annotation
                sample_run_to_sample_name = {}
                with open(sample_annotation_path, 'r') as f:
                    reader = csv.DictReader(f)
                    for row in reader:
                        sample_name = row.get('sample_name', '').strip()
                        run = row.get('run', '').strip()
                        if sample_name and run:
                            sample_run = f"{sample_name}_{run}"
                            sample_run_to_sample_name[sample_run] = sample_name
            except Exception as e:
                log.warning('Could not load annotation file for prealign stats mapping: {}'.format(str(e)))
                sample_run_to_sample_name = {}
            
            # Parse each prealign stats file
            for f in prealign_stats_files:
                # Extract sample_run from filename (e.g., "test_1.prealign.stats.tsv" -> "test_1")
                sample_run = f['s_name'].replace('.prealign.stats', '')
                # Map to sample_name
                sample_name = sample_run_to_sample_name.get(sample_run, sample_run.split('_')[0] if '_' in sample_run else sample_run)
                
                # Parse the stats file to extract percent_filtered from "total" row
                prealign_data = self.parse_prealign_stats(f['f'])
                if prealign_data and 'percent_filtered' in prealign_data:
                    if sample_name not in self.atacseq_data:
                        self.atacseq_data[sample_name] = {}
                    # Store percent_filtered (will aggregate if multiple runs per sample)
                    if 'percent_filtered' not in self.atacseq_data[sample_name]:
                        self.atacseq_data[sample_name]['percent_filtered'] = []
                    self.atacseq_data[sample_name]['percent_filtered'].append(prealign_data['percent_filtered'])
            
            # Aggregate percent_filtered per sample (average if multiple runs)
            for sample_name in self.atacseq_data:
                if 'percent_filtered' in self.atacseq_data[sample_name] and isinstance(self.atacseq_data[sample_name]['percent_filtered'], list):
                    values = self.atacseq_data[sample_name]['percent_filtered']
                    if values:
                        avg_percent = sum(values) / len(values)
                        self.atacseq_data[sample_name]['percent_filtered'] = avg_percent
                    else:
                        del self.atacseq_data[sample_name]['percent_filtered']
            
            log.info('Found prealign stats file for {} sample_runs'.format(len(prealign_stats_files)))
        
        # Remove ignored samples if there is any
        self.atacseq_tss_data = self.ignore_samples(self.atacseq_tss_data)
        # Remove ignored samples if there is any
        self.atacseq_data = self.ignore_samples(self.atacseq_data)

        # Load the sample annotation sheet
        sample_sas_path = config.annotation
        sample_sas = csv.DictReader(open(sample_sas_path, 'r'))
        self.sample_sas_dict = {}
        for k in sample_sas:
            print(k)
            self.sample_sas_dict[k['sample_name']] = k

        # Check if there are any paired end sample in the current project
        self.pairedSampleExists = False
        for sample in self.sample_sas_dict:
            if self.sample_sas_dict[sample]['read_type'] == 'paired':
                self.pairedSampleExists = True
        
        # Get the genome version
        self.genome_version = config.genome

        # Parse reproducibility and BAM correlation stats
        self.reproducibility_stats = self.parse_reproducibility_stats()
        self.bam_correlation_stats = self.parse_bam_correlation_stats()
        
        # Map replicate groups to samples for general stats
        self.replicate_group_to_samples = self.build_replicate_group_mapping()

        # Add stats to general table
        self.add_atacseq_to_general_stats()

        # Add download links table
        self.add_download_table()

        # Add TSS line graph
        self.add_tss_plot()

    def parse_atacseq_stats(self, f):
        data = {}
        for l in f.splitlines():
            s = l.split('\t')
            data[s[0]] = s[1]
        return data

    def parse_atacseq_tss(self, f):
        data = OrderedDict()
        count = 0
        max_value = 0.0
        for l in f.splitlines():
            s = l.split(',')
            if s[0] == 'base':
                continue
            if float(s[1]) > max_value:
                max_value = float(s[1])
            count += 1
            if count % 10 == 0:
                data[int(s[0])] = float(s[1])
        return data, max_value

    def parse_align_stats(self, f):
        """
        Parse align_stats file which contains mitochondrial_fraction.
        Format: metric\tvalue (single line)
        """
        data = {}
        for l in f.splitlines():
            l = l.strip()
            if not l:  # Skip empty lines
                continue
            s = l.split('\t')
            if len(s) >= 2:
                metric = s[0].strip()
                value = s[1].strip()
                data[metric] = value
        return data

    def parse_prealign_stats(self, f):
        """
        Parse prealign stats file to extract percent_filtered from the "total" row.
        Format: prealignment\treads_before\treads_after\treads_filtered\tpercent_filtered
        """
        data = {}
        lines = f.splitlines()
        if len(lines) < 2:
            return data
        
        # Skip header line
        for line in lines[1:]:
            line = line.strip()
            if not line:
                continue
            fields = line.split('\t')
            if len(fields) >= 5:
                prealignment = fields[0].strip()
                if prealignment == 'total':
                    try:
                        percent_filtered = float(fields[4].strip())
                        data['percent_filtered'] = percent_filtered
                    except (ValueError, IndexError):
                        pass
                    break
        return data

    def parse_reproducibility_stats(self):
        """
        Parse reproducibility QC stats TSV file.
        Format: Sample\tReproducibility\tN_Optimal_Peaks\tRescue_Ratio\tSelf_Consistency_Ratio
        """
        stats = {}
        try:
            import glob
            # Search in common report directory locations
            search_paths = []
            
            # Try config.data_dir first
            if hasattr(config, 'data_dir') and config.data_dir:
                search_paths.append(config.data_dir)
            
            # Try report directory relative to data_dir
            if hasattr(config, 'data_dir') and config.data_dir:
                report_dir = os.path.join(config.data_dir, 'report')
                if os.path.exists(report_dir):
                    search_paths.append(report_dir)
            
            # Try current working directory
            search_paths.append(os.getcwd())
            
            # Search for the file
            for base_dir in search_paths:
                pattern = os.path.join(base_dir, '**', 'reproducibility_stats_mqc.tsv')
                files = glob.glob(pattern, recursive=True)
                if files:
                    filepath = files[0]  # Use first match
                    log.info('Found reproducibility stats file: {}'.format(filepath))
                    with open(filepath, 'r') as f:
                        lines = f.readlines()
                        if len(lines) < 2:
                            continue
                        # Skip header
                        for line in lines[1:]:
                            parts = line.strip().split('\t')
                            if len(parts) >= 5:
                                group = parts[0]
                                reproducibility = parts[1]
                                n_optimal = parts[2] if parts[2] != 'N/A' else None
                                rescue_ratio = parts[3] if parts[3] != 'N/A' else None
                                self_consistency_ratio = parts[4] if parts[4] != 'N/A' else None
                                
                                stats[group] = {
                                    'reproducibility': reproducibility,
                                    'n_optimal_peaks': int(n_optimal) if n_optimal and n_optimal.replace('.', '').isdigit() else None,
                                    'rescue_ratio': float(rescue_ratio) if rescue_ratio and rescue_ratio != 'N/A' else None,
                                    'self_consistency_ratio': float(self_consistency_ratio) if self_consistency_ratio and self_consistency_ratio != 'N/A' else None,
                                }
                    break
        except Exception as e:
            log.warning('Could not parse reproducibility stats: {}'.format(str(e)))
        return stats

    def parse_bam_correlation_stats(self):
        """
        Parse BAM correlation stats TSV file.
        Format: Sample\tBAM_Correlation\tCorrelation_Quality\tN_Replicates
        """
        stats = {}
        try:
            import glob
            # Search in common report directory locations
            search_paths = []
            
            # Try config.data_dir first
            if hasattr(config, 'data_dir') and config.data_dir:
                search_paths.append(config.data_dir)
            
            # Try report directory relative to data_dir
            if hasattr(config, 'data_dir') and config.data_dir:
                report_dir = os.path.join(config.data_dir, 'report')
                if os.path.exists(report_dir):
                    search_paths.append(report_dir)
            
            # Try current working directory
            search_paths.append(os.getcwd())
            
            # Search for the file
            for base_dir in search_paths:
                pattern = os.path.join(base_dir, '**', 'bam_correlation_stats_mqc.tsv')
                files = glob.glob(pattern, recursive=True)
                if files:
                    filepath = files[0]  # Use first match
                    log.info('Found BAM correlation stats file: {}'.format(filepath))
                    with open(filepath, 'r') as f:
                        lines = f.readlines()
                        if len(lines) < 2:
                            continue
                        # Skip header
                        for line in lines[1:]:
                            parts = line.strip().split('\t')
                            if len(parts) >= 4:
                                group = parts[0]
                                bam_correlation = parts[1] if parts[1] != 'N/A' else None
                                quality = parts[2]
                                n_replicates = parts[3] if parts[3].isdigit() else None
                                
                                stats[group] = {
                                    'bam_correlation': float(bam_correlation) if bam_correlation and bam_correlation != 'N/A' else None,
                                    'quality': quality,
                                    'n_replicates': int(n_replicates) if n_replicates else None,
                                }
                    break
        except Exception as e:
            log.warning('Could not parse BAM correlation stats: {}'.format(str(e)))
        return stats

    def build_replicate_group_mapping(self):
        """
        Build a mapping from replicate group names to individual sample names.
        Uses the annotation file to determine which samples belong to which replicate group.
        """
        mapping = {}
        try:
            sample_annotation_path = config.annotation
            with open(sample_annotation_path, 'r') as f:
                reader = csv.DictReader(f)
                for row in reader:
                    sample_name = row.get('sample_name', '').strip()
                    replicate_group = row.get('replicate_sample_name', '').strip()
                    if sample_name and replicate_group:
                        if replicate_group not in mapping:
                            mapping[replicate_group] = []
                        mapping[replicate_group].append(sample_name)
        except Exception as e:
            log.warning('Could not build replicate group mapping: {}'.format(str(e)))
        return mapping

    def add_atacseq_to_general_stats(self):
        data = {}
        for sample_name in self.atacseq_data:
            data[sample_name] = {}
            if hasattr(config, 'exploratory_columns'):
                for column in config.exploratory_columns:
                    if column in self.sample_sas_dict[sample_name]:
                        data[sample_name][column] = self.sample_sas_dict[sample_name][column]
            if 'NSC' in self.atacseq_data[sample_name] and self.atacseq_data[sample_name]['NSC'] != 'nan':
                try:
                    value = float(self.atacseq_data[sample_name]['NSC'])
                except ValueError as err:
                    print(err)
                    value = 'NaN'
                data[sample_name]['NSC'] = value
            if 'RSC' in self.atacseq_data[sample_name] and self.atacseq_data[sample_name]['RSC'] != 'nan':
                try:
                    value = float(self.atacseq_data[sample_name]['RSC'])
                except:
                    value = 'NaN'
                data[sample_name]['RSC'] = value
            if 'peaks' in self.atacseq_data[sample_name]:
                try:
                    value = int(self.atacseq_data[sample_name]['peaks'])
                except:
                    value = None
                data[sample_name]['peaks'] = value
            if 'filtered_peaks' in self.atacseq_data[sample_name]:
                try:
                    value = int(self.atacseq_data[sample_name]['filtered_peaks'])
                except:
                    value = None
                data[sample_name]['filtered_peaks'] = value
            if 'frip' in self.atacseq_data[sample_name]:
                try:
                    value = float(self.atacseq_data[sample_name]['frip'])
                except:
                    value = None
                data[sample_name]['frip'] = value
            if 'regulatory_fraction' in self.atacseq_data[sample_name]:
                try:
                    value = float(self.atacseq_data[sample_name]['regulatory_fraction'])
                except:
                    value = None
                data[sample_name]['regulatory_fraction'] = value
            if 'tss_max' in self.atacseq_data[sample_name]:
                data[sample_name]['tss_max'] = self.atacseq_data[sample_name]['tss_max']
            if 'mitochondrial_fraction' in self.atacseq_data[sample_name]:
                try:
                    value = float(self.atacseq_data[sample_name]['mitochondrial_fraction'])
                    data[sample_name]['mitochondrial_fraction'] = value
                except (ValueError, TypeError):
                    data[sample_name]['mitochondrial_fraction'] = None
            
            if 'n_tot' in self.atacseq_data[sample_name]:
                try:
                    value = int(self.atacseq_data[sample_name]['n_tot'])
                except:
                    value = None
                data[sample_name]['n_tot'] = value
            
            if 'n_nondups' in self.atacseq_data[sample_name]:
                try:
                    value = int(self.atacseq_data[sample_name]['n_nondups'])
                except:
                    value = None
                data[sample_name]['n_nondups'] = value
            
            if 'nrf' in self.atacseq_data[sample_name]:
                try:
                    value = float(self.atacseq_data[sample_name]['nrf'])
                except:
                    value = None
                data[sample_name]['nrf'] = value
            
            if 'percent_filtered' in self.atacseq_data[sample_name]:
                try:
                    value = float(self.atacseq_data[sample_name]['percent_filtered'])
                except:
                    value = None
                data[sample_name]['percent_filtered'] = value
            
            # Add reproducibility metrics if this sample belongs to a replicate group
            for group_name, samples_in_group in self.replicate_group_to_samples.items():
                if sample_name in samples_in_group:
                    # Add reproducibility stats
                    if group_name in self.reproducibility_stats:
                        rep_stats = self.reproducibility_stats[group_name]
                        if rep_stats.get('rescue_ratio') is not None:
                            data[sample_name]['rescue_ratio'] = rep_stats['rescue_ratio']
                        if rep_stats.get('self_consistency_ratio') is not None:
                            data[sample_name]['self_consistency_ratio'] = rep_stats['self_consistency_ratio']
                        if rep_stats.get('n_optimal_peaks') is not None:
                            data[sample_name]['n_optimal_peaks'] = rep_stats['n_optimal_peaks']
                        data[sample_name]['reproducibility_status'] = rep_stats.get('reproducibility', 'unknown')
                    
                    # Add BAM correlation stats
                    if group_name in self.bam_correlation_stats:
                        bam_stats = self.bam_correlation_stats[group_name]
                        if bam_stats.get('bam_correlation') is not None:
                            data[sample_name]['bam_correlation'] = bam_stats['bam_correlation']
                    break  # Sample can only belong to one group
        
        headers = OrderedDict()
        if hasattr(config, 'exploratory_columns'):
            for column in config.exploratory_columns:
                log.info('Adding exploratory column {}'.format(column))
                headers[column] = {
                    'description': column,
                    'title': column,
                    'scale': False }
        else:
            log.warning("No exploratory columns were specified in the config")
        
        # Get list of columns to show/hide from config (if provided)
        visible_columns = getattr(config, 'atacseq_general_stats_columns', None)
        if visible_columns is not None:
            log.info('Using config-specified columns for General Statistics: {}'.format(visible_columns))
        
        # Helper function to conditionally add columns
        def add_header_if_visible(key, header_config):
            """Add header only if visible_columns is None (all visible) or key is in visible_columns"""
            if visible_columns is None or key in visible_columns:
                headers[key] = header_config
            else:
                # Still add but mark as hidden
                header_config['hidden'] = True
                headers[key] = header_config

        add_header_if_visible('peaks', {
            'description': 'Number of detected peaks',
            'title': 'Peaks',
            'scale': 'Greens',
            'format': '{:,.0f}'
        })
        
        add_header_if_visible('filtered_peaks', {
            'description': 'Number of peaks remaining after filtering',
            'title': 'Filtered\nPeaks',
            'scale': 'Greens',
            'format': '{:,.0f}'
        })
        
        add_header_if_visible('NSC', {
            'description': 'Normalized Strand Cross-correlation Coefficient',
            'title': 'NSC',
            'scale': 'Reds',
            'min': 0.0,
            'max': 2.0,
            'format': '{:.2f}'
        })
        
        add_header_if_visible('NSC_PCT', {
            'description': 'NSC Percentile Among All ATAC-seq samples',
            'title': 'NSC_PCT',
            'scale': 'Reds',
            'suffix': '%',
            'max': 100,
            'format': '{:,.0f}'
        })
        
        add_header_if_visible('RSC', {
            'description': 'Relative Strand Cross-correlation Coefficient',
            'title': 'RSC',
            'scale': 'Reds',
            'min': 0.0,
            'max': 2.0,
            'format': '{:.2f}'
        })
        
        add_header_if_visible('RSC_PCT', {
            'description': 'RSC Percentile Among All ATAC-seq samples',
            'title': 'RSC_PCT',
            'scale': 'Reds',
            'suffix': '%',
            'max': 100,
            'format': '{:,.0f}'
        })
        
        add_header_if_visible('frip', {
            'description': 'Fraction of Reads in Peaks',
            'title': 'FRiP',
            'scale': 'Reds-rev',
            'min': 0.0,
            'max': 1.0,
            'format': '{:.2f}'
        })
        
        add_header_if_visible('regulatory_fraction', {
            'description': 'Fraction of Reads in Regulatory Regions',
            'title': 'Regulatory',
            'scale': 'Reds-rev',
            'min': 0.0,
            'max': 1.0,
            'format': '{:.2f}'
        })
        
        add_header_if_visible('tss_max', {
            'description': 'The peak value of TSS enrichment',
            'title': 'TSS',
            'scale': 'Reds-rev',
            'format': '{:.1f}'
        })
        
        add_header_if_visible('mitochondrial_fraction', {
            'description': 'Fraction of Reads from Mitochondria',
            'title': 'Mitochondrial DNA',
            'scale': 'Reds-rev',
            'min': 0.0,
            'max': 1.0,
            'format': '{:.2f}'
        })
        
        add_header_if_visible('n_tot', {
            'description': 'Total number of alignments processed by Samblaster',
            'title': 'Total\nAlignments',
            'scale': 'Blues',
            'format': '{:,.0f}'
        })
        
        add_header_if_visible('n_nondups', {
            'description': 'Number of non-duplicate alignments from Samblaster',
            'title': 'Non-Duplicate\nAlignments',
            'scale': 'Greens',
            'format': '{:,.0f}'
        })
        
        add_header_if_visible('nrf', {
            'description': 'Non-Redundant Fraction (1 - duplication rate)',
            'title': 'NRF',
            'scale': 'Greens',
            'min': 0.0,
            'max': 1.0,
            'format': '{:.4f}'
        })
        
        add_header_if_visible('percent_filtered', {
            'description': 'Percentage of reads filtered during prealignment',
            'title': 'Prealign\nFiltered',
            'scale': 'Reds-rev',
            'min': 0.0,
            'max': 1.0,
            'format': '{:.2f}'
        })
        
        # Reproducibility metrics (for replicate groups)
        add_header_if_visible('rescue_ratio', {
            'description': 'Rescue Ratio: max(Np, Nt) / min(Np, Nt). Should be ≤2.0 (ENCODE standard)',
            'title': 'Rescue\nRatio',
            'scale': 'RdYlGn',
            'min': 1.0,
            'max': 5.0,
            'format': '{:.3f}'
        })
        
        add_header_if_visible('self_consistency_ratio', {
            'description': 'Self-Consistency Ratio: max(N1, N2) / min(N1, N2). Should be ≤2.0 (ENCODE standard)',
            'title': 'Self-Consist.\nRatio',
            'scale': 'RdYlGn',
            'min': 1.0,
            'max': 5.0,
            'format': '{:.3f}'
        })
        
        add_header_if_visible('n_optimal_peaks', {
            'description': 'Number of peaks in optimal set (selected based on reproducibility)',
            'title': 'Optimal\nPeaks',
            'scale': 'Greens',
            'format': '{:,.0f}'
        })
        
        add_header_if_visible('reproducibility_status', {
            'description': 'Reproducibility status: pass/borderline/fail (ENCODE standard)',
            'title': 'Reproducibility',
            'scale': False,
            'format': '{:s}'
        })
        
        # BAM correlation metrics (for replicate groups)
        add_header_if_visible('bam_correlation', {
            'description': 'Mean Spearman correlation between BAM files (off-diagonal mean). Higher is better.',
            'title': 'BAM\nCorrelation',
            'scale': 'RdYlGn',
            'min': 0.0,
            'max': 1.0,
            'format': '{:.4f}'
        })
        
        self.general_stats_addcols(data, headers)

    def add_download_table(self):
        # Create a table with download links to various files
        results_url = '../results' #os.path.join(config.base_url, config.project_uuid, 'results')
#         project_url = os.path.join(config.base_url, config.project_uuid)

        # Configuration for the MultiQC table
        table_config = {
            'namespace': 'Download links',  # Name for grouping. Prepends desc and is in Config Columns modal
            'id': 'download_links',  # ID used for the table
            'table_title': 'Download ATAC-seq data',  # Title of the table. Used in the column config modal
            'save_file': False,  # Whether to save the table data to a file
            'sortRows': False,  # Whether to sort rows alphabetically
            'col1_header': 'Sample Name',  # The header used for the first column
            'no_beeswarm': True,
            'scale': False
        }

        # Configuration for the header row
        headers = OrderedDict()
        headers['BAM'] = {
            'title': 'BAM',
            'description': 'Bowtie2 alignment results in BAM format.',
            'scale': False,
            'hidden': False
        }
        headers['filtered_BAM'] = {
            'title': 'Filtered BAM',
            'description': 'BAM files without the low quality alignments.',
            'scale': False,
            'hidden': False
        }
        headers['filtered_peaks'] = {
            'title': 'Peaks',
            'description': 'Peak calls by MACS2 in bed format.',
            'scale': False,
            'hidden': False
        }
        headers['summits_bed'] = {
            'title': 'Summits',
            'description': 'Summits of the peaks in bed format.',
            'scale': False,
            'hidden': False
        }
        headers['motifs'] = {
            'title': 'Motifs',
            'description': 'HOMER motif analysis results.',
            'scale': False,
            'hidden': False
        }
#         headers['coverage_bigwig'] = {
#             'title': 'Coverage BigWig',
#             'description': 'Genome wide coverage data in UCSC bigWig format.',
#             'scale': False,
#             'hidden': False
#         }

        # Fill the download table with URLs
#         igv_links = []
        sample_names = []
        data = OrderedDict()
        for sample_name in self.atacseq_data:
            sample_names.append(sample_name)
            # generate links list for loading them on IGV
#             igv_link = project_url + '/hub/' + self.genome_version + '/' + sample_name + '.bigWig'
#             igv_links.append(igv_link)
            sample_bam_url = '{}/{}/mapped/{}.bam'.format(results_url, sample_name, sample_name)
            sample_bai_url = '{}/{}/mapped/{}.bam.bai'.format(results_url, sample_name, sample_name)
            sample_filtered_bam_url = '{}/{}/mapped/{}.filtered.bam'.format(results_url, sample_name, sample_name)
            sample_filtered_bai_url = '{}/{}/mapped/{}.filtered.bam.bai'.format(results_url, sample_name, sample_name)
            sample_peaks_url = '{}/{}/peaks/{}_peaks.narrowPeak'.format(results_url, sample_name, sample_name)
            sample_annotated_peaks_url = '{}/{}/peaks/{}_peaks.narrowPeak.annotated.tsv'.format(results_url, sample_name, sample_name)
            sample_summits_url = '{}/{}/peaks/{}_summits.bed'.format(results_url, sample_name, sample_name)
            sample_known_motifs_url = '{}/{}/homer/knownResults.html'.format(results_url, sample_name)
            sample_denovo_motifs_url = '{}/{}/homer/homerResults.html'.format(results_url, sample_name)
#             sample_bigwig_url = '../hub/{}.bigWig'.format(sample_name)
            data[sample_name] = {
                'BAM': '<a href={}>{} BAM</a></br><a href={}>{} BAI</a>'.format(sample_bam_url, sample_name, sample_bai_url, sample_name),
                'filtered_BAM': '<a href={}>{} flt BAM</a></br><a href={}>{} flt BAI</a>'.format(sample_filtered_bam_url, sample_name, sample_filtered_bai_url, sample_name),
                'filtered_peaks': '<a href={}>{} Peaks</a></br><a href={}>{} Annotated Peaks</a>'.format(sample_peaks_url, sample_name, sample_annotated_peaks_url, sample_name),
                'summits_bed': '<a href={}>{} Summits</a>'.format(sample_summits_url, sample_name),
#                 'coverage_bigwig': '<a href={}>{} bigWig</a>'.format(sample_bigwig_url, sample_name),
                'motifs': '<a href={}>{} Known</a></br><a href={}>{} DeNovo</a>'.format(sample_known_motifs_url, sample_name, sample_denovo_motifs_url, sample_name)
            }

#         # Generate the UCSC genome browser link
#         track_hubs_url = project_url + '/hub/hub.txt'
#         genome_browser_url = 'http://genome-euro.ucsc.edu/cgi-bin/hgTracks?db={}&hubUrl={}'.format(self.genome_version, track_hubs_url)
#         section_description = '<a href={} target=_blank>Click here to view the coverage tracks on UCSC Genome Browser <span class="glyphicon glyphicon-new-window"></span></a>'.format(genome_browser_url)

#         igv_load_link = '<a href=http://localhost:60151/load?file={}?names={}?genome={}?goto=chr1>Click here to load the coverage tracks on your IGV session (Genome: {}) <span class="glyphicon glyphicon-new-window"></span></a>'.format(
#             ','.join(igv_links),
#             ','.join(sample_names),
#             self.genome_version,
#             self.genome_version
#         )
#         section_description += '<br>' + igv_load_link
        # Finally add a MultiQC section together with the URL table
        self.add_section(
            name='Download Links & Coverage Tracks',
            anchor='atacseq_download',
#             description=section_description,
            helptext='You can click on the table elements to download the files.',
            plot=table.plot(data, headers, table_config)
        )

    def add_tss_plot(self):
        tss_plot_config = {
            'xlab': 'Distance from TSS in Bp',
            'ylab': 'Normalized Coverage'
        }

        self.add_section(
            name='Coverage around TSS',
            anchor='atacseq_tss',
            description='Coverage plot of sites around TSS sites',
            helptext='This plot shows the aggregated and normalized coverage around the transcription start sites (TSS)',
            plot=linegraph.plot(data=self.atacseq_tss_data, pconfig=tss_plot_config)
        )
