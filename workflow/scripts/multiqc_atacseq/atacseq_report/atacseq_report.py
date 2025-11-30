#!/usr/bin/env python
""" MultiQC BSF plugin functions
We can add any custom Python functions here and call them
using the setuptools plugin hooks.
"""


from __future__ import print_function
from pkg_resources import get_distribution
import logging
import os, csv

from multiqc.utils import report, util_functions, config

# Initialise the main MultiQC logger
log = logging.getLogger('multiqc')

# Save this plugin's version number (defined in setup.py) to the MultiQC config
config.atacseq_report_version = get_distribution("atacseq_report").version


# Add default config options for the things that are used in atacseq_report
def atacseq_report_execution_start():
    """
    Code to execute after the config files and
    command line flags have been parsed self.
    this setuptools hook is the earliest that will be able
    to use custom command line flags.
    """
    # Halt execution if we've disabled the plugin
    if config.kwargs.get('disable_atacseq_report', True):
        return None

    log.info("Running atacseq_report MultiQC Plugin v{}, use --disable-atacseq-report to disable".format(config.atacseq_report_version))

    # Add to the search patterns used by atacseq module
    if 'atacseq' not in config.sp:
        config.update_dict(config.sp, {'atacseq': {'fn': '*.stats.tsv', 'contents': 'frip'}})
        log.info("updated config.sp for atacseq")
    if 'atacseq/tss' not in config.sp:
        config.update_dict(config.sp, {'atacseq/tss': {'fn': '*TSS.csv', 'contents': 'count'}})
    if 'atacseq/align_stats' not in config.sp:
        config.update_dict(config.sp, {'atacseq/align_stats': {'fn': '*.align.stats.tsv', 'contents': 'mitochondrial_fraction'}})
        log.info("updated config.sp for atacseq/align_stats")
    else:
        log.info("atacseq/align_stats search pattern already exists in config.sp")


def atacseq_report_after_modules():
    """
    Hook that runs after all modules are initialized.
    This allows us to access data from other modules like Sambamba/Samblaster.
    """
    # This hook runs after modules, but we'll handle data access in the module itself
    # using a different approach
    pass
