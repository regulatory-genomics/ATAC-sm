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
                                            href='https://github.com/epigen/atacseq_pipeline',
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
        log.info('Searching for align_stats files with pattern atacseq/align_stats')
        # Check if search pattern exists
        if 'atacseq/align_stats' not in config.sp:
            log.warning('Search pattern atacseq/align_stats not found in config.sp! Available patterns: {}'.format(list(config.sp.keys())))
        else:
            log.debug('Search pattern atacseq/align_stats found in config.sp: {}'.format(config.sp.get('atacseq/align_stats')))
        
        align_stats_files = list(self.find_log_files(sp_key='atacseq/align_stats'))
        log.info('Found {} align_stats files using search pattern'.format(len(align_stats_files)))
        
        if len(align_stats_files) == 0:
            # Try to find files manually for debugging
            log.warning('No align_stats files found with search pattern. Checking all found files...')
            all_files = list(self.find_log_files())
            align_files = [f for f in all_files if 'align.stats' in f.get('fn', '')]
            log.warning('Found {} files with "align.stats" in filename: {}'.format(
                len(align_files), [f.get('fn', '') for f in align_files]))
            # Also try without the search pattern key
            log.warning('Trying to find files with pattern *.align.stats.tsv directly...')
            direct_files = [f for f in all_files if f.get('fn', '').endswith('.align.stats.tsv')]
            log.warning('Found {} files ending with .align.stats.tsv: {}'.format(
                len(direct_files), [f.get('fn', '') for f in direct_files]))
        
        for f in align_stats_files:
            # Extract sample name from filename (e.g., "sample.align.stats.tsv" -> "sample")
            # The filename is like "test.align.stats.tsv", extract "test"
            s_name = f['s_name']
            if s_name.endswith('.align.stats'):
                sample_name = s_name.replace('.align.stats', '')
            elif s_name.endswith('.align'):
                sample_name = s_name.replace('.align', '')
            else:
                sample_name = s_name
            
            log.debug('Processing align_stats file: {} -> s_name: {} -> sample: {}'.format(f.get('fn', ''), s_name, sample_name))
            align_stats_data = self.parse_align_stats(f['f'])
            # Merge align_stats data into atacseq_data
            if sample_name not in self.atacseq_data:
                self.atacseq_data[sample_name] = {}
            self.atacseq_data[sample_name].update(align_stats_data)
            log.debug('Parsed align_stats for {}: {}'.format(sample_name, align_stats_data))
        log.info('Found align_stats file for {} ATAC-seq samples'.format(len(align_stats_files)))
        
        # Parse Samblaster logs for n_tot and n_nondups
        # Samblaster logs contain: "Marked X of Y (Z%) read ids as duplicates"
        # n_dups = X, n_tot = Y, n_nondups = Y - X
        log.info('Searching for Samblaster log files...')
        samblaster_files = list(self.find_log_files('samblaster'))
        log.info('Found {} Samblaster log files'.format(len(samblaster_files)))
        
        for f in samblaster_files:
            sample_name = f['s_name']
            log.debug('Processing Samblaster log for sample: {}'.format(sample_name))
            
            # Parse the log file
            for line in f['f'].splitlines():
                if 'Marked' in line and 'read ids as duplicates' in line:
                    # Example: "samblaster: Marked 3894 of 7458 (52.21%) read ids as duplicates using 2864k memory in 0.014S CPU seconds and 22S wall time."
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
                        
                        log.debug('Parsed Samblaster for {}: n_tot={}, n_nondups={}, n_dups={}, nrf={:.4f}'.format(
                            sample_name, n_tot, n_nondups, n_dups, self.atacseq_data[sample_name]['nrf']))
                        break
        
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

        # Add stats to general table
        # Note: This is called during module initialization, so other modules' data
        # might not be available yet. We'll try to access it, but if not available,
        # the columns will still be added (they just won't have data)
        log.info('Calling add_atacseq_to_general_stats()...')
        self.add_atacseq_to_general_stats()
        log.info('Finished add_atacseq_to_general_stats()')

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

    def add_atacseq_to_general_stats(self):
        # Get Sambamba/Samblaster data from MultiQC's built-in module (n_tot, n_nondups)
        # MultiQC automatically parses these logs - access via the general_stats_data dict
        # which is populated by BaseMultiqcModule
        sambamba_data = None
        module_name = None
        
        # Access general_stats_data from the report object
        # In MultiQC, this is populated as modules add their data
        try:
            from multiqc.utils import report
            log.debug('Checking report.general_stats_data availability...')
            # Try accessing general_stats_data directly
            if hasattr(report, 'general_stats_data'):
                log.debug('report.general_stats_data exists: {}'.format(bool(report.general_stats_data)))
                if report.general_stats_data:
                    available_modules = list(report.general_stats_data.keys())
                    log.info('Available modules in general_stats_data: {}'.format(available_modules))
                    # Check for sambamba module (Sambamba markdup)
                    if 'sambamba' in report.general_stats_data:
                        sambamba_data = report.general_stats_data['sambamba']
                        module_name = 'sambamba'
                        log.info('Found Sambamba module data with {} samples'.format(len(sambamba_data)))
                    # Check for samblaster module (alternative)
                    elif 'samblaster' in report.general_stats_data:
                        sambamba_data = report.general_stats_data['samblaster']
                        module_name = 'samblaster'
                        log.info('Found Samblaster module data with {} samples'.format(len(sambamba_data)))
                        # Debug: show sample names and available keys
                        for sample in list(sambamba_data.keys())[:3]:  # Show first 3 samples
                            log.debug('Samblaster sample {} has keys: {}'.format(sample, list(sambamba_data[sample].keys())))
                    else:
                        log.warning('Neither sambamba nor samblaster found in general_stats_data')
                else:
                    log.warning('report.general_stats_data is empty or None')
            else:
                log.warning('report.general_stats_data attribute does not exist')
        except Exception as e:
            log.warning('Could not access report.general_stats_data: {}'.format(str(e)))
            import traceback
            log.debug(traceback.format_exc())
        
        # Extract n_tot and n_nondups from Sambamba/Samblaster data
        if sambamba_data:
            for sample_name in sambamba_data:
                if sample_name not in self.atacseq_data:
                    self.atacseq_data[sample_name] = {}
                sample_data = sambamba_data[sample_name]
                # Extract n_tot and n_nondups
                if 'n_tot' in sample_data:
                    self.atacseq_data[sample_name]['n_tot'] = sample_data['n_tot']
                    log.debug('Added n_tot={} for sample {}'.format(sample_data['n_tot'], sample_name))
                if 'n_nondups' in sample_data:
                    self.atacseq_data[sample_name]['n_nondups'] = sample_data['n_nondups']
                    log.debug('Added n_nondups={} for sample {}'.format(sample_data['n_nondups'], sample_name))
            num_samples = len([s for s in sambamba_data if 'n_tot' in sambamba_data.get(s, {})])
            log.info('Extracted n_tot/n_nondups for {} samples from {} module'.format(num_samples, module_name))
        else:
            log.warning('Sambamba/Samblaster module data not found - n_tot/n_nondups columns will be empty')
        
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
                    log.debug('Added mitochondrial_fraction={} for sample {}'.format(value, sample_name))
                except (ValueError, TypeError) as e:
                    log.warning('Could not convert mitochondrial_fraction to float for sample {}: {}'.format(sample_name, str(e)))
                    data[sample_name]['mitochondrial_fraction'] = None
            else:
                log.debug('mitochondrial_fraction not found for sample: {} (available keys: {})'.format(
                    sample_name, list(self.atacseq_data[sample_name].keys())))
            
            if 'n_tot' in self.atacseq_data[sample_name]:
                try:
                    value = int(self.atacseq_data[sample_name]['n_tot'])
                except:
                    value = None
                data[sample_name]['n_tot'] = value
            else:
                log.debug('n_tot not found for sample: {}'.format(sample_name))
            
            if 'n_nondups' in self.atacseq_data[sample_name]:
                try:
                    value = int(self.atacseq_data[sample_name]['n_nondups'])
                except:
                    value = None
                data[sample_name]['n_nondups'] = value
            else:
                log.debug('n_nondups not found for sample: {}'.format(sample_name))
            
            if 'nrf' in self.atacseq_data[sample_name]:
                try:
                    value = float(self.atacseq_data[sample_name]['nrf'])
                except:
                    value = None
                data[sample_name]['nrf'] = value
            else:
                log.debug('nrf not found for sample: {}'.format(sample_name))
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
            'description': 'Total number of alignments processed by Sambamba markdup',
            'title': 'Total\nAlignments',
            'scale': 'Blues',
            'format': '{:,.0f}'
        })
        
        add_header_if_visible('n_nondups', {
            'description': 'Number of non-duplicate alignments from Sambamba markdup',
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
