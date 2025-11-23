#!/bin/env python

#### libraries
import os
import pandas as pd
import pybedtools as bedtools

#### configurations

# input
peakfiles = snakemake.input["summits_bed"]
blacklist_file = snakemake.input["blacklisted_regions"]
chrom_file = snakemake.input["chromosome_sizes"]

# output
consensus_regions_path = snakemake.output["consensus_regions"]

# parameters
slop_extension = snakemake.params["slop_extension"]

# load summits and generate consensus regions using (py)bedtools
output_bed = None

for peakfile in peakfiles:
    # Check if file exists and is not empty
    if os.path.exists(peakfile) and os.path.getsize(peakfile) > 0:
        peak_bed = bedtools.BedTool(peakfile)
        if (blacklist_file is not None):
            peak_bed=peak_bed.intersect(blacklist_file,v=True, wa=True)

        peak_bed = peak_bed.slop(g=chrom_file, b=slop_extension)

        if (output_bed is None):
            output_bed = peak_bed
        else:
            output_bed = output_bed.cat(peak_bed,force_truncate=True)

# Handle empty output
if output_bed is None or len(output_bed) == 0:
    # Create empty BED file
    with open(consensus_regions_path, 'w') as f:
        pass  # Create empty file
else:
    output_bed.saveas(consensus_regions_path)
    peaks = bedtools.BedTool(consensus_regions_path).sort(faidx=chrom_file).to_dataframe(names=['CHR','START','END'],dtype={'START':int,'END':int})
    
    # Create ID column with proper formatting
    peaks['ID'] = peaks.index.to_series().apply(lambda x: "CONS{:011d}".format(x))
    
    # Reorder columns to match BED format: chrom, start, end, name
    peaks_bed = peaks[['CHR','START','END','ID']].copy()
    peaks_bed.columns = ['chrom','start','end','name']
    bedtools.BedTool().from_dataframe(peaks_bed).saveas(consensus_regions_path)
