#!/bin/env python

#### libraries
import os
import pandas as pd
import pybedtools as bedtools

#### configurations

# input
consensus_regions_path = snakemake.input["consensus_regions"]
peakfile_path = snakemake.input["peakfile"]
chrom_file = snakemake.input["chromosome_sizes"]

# output
quant_support_path = snakemake.output["quant_support"]

# parameters
sample = snakemake.wildcards["sample"]

# Check if consensus_regions file exists and is not empty
if not os.path.exists(consensus_regions_path) or os.path.getsize(consensus_regions_path) == 0:
    # Create empty result file
    result = pd.DataFrame(columns=[sample])
    result.T.to_csv(quant_support_path)
else:
    # quantify peak support within sample
    consensus_peaks = bedtools.BedTool(consensus_regions_path)
    consensus_peaks_df = consensus_peaks.to_dataframe()
    
    # Check if 'name' column exists, if not create it from index
    if 'name' not in consensus_peaks_df.columns:
        if len(consensus_peaks_df) > 0:
            # Create name column from index if it doesn't exist
            consensus_peaks_df['name'] = consensus_peaks_df.index.to_series().apply(lambda x: f"CONS{x:011d}")
        else:
            consensus_peaks_df['name'] = []
    
    consensus_peaks_df = consensus_peaks_df.set_index('name')
    
    # Initialize result dataframe
    result = pd.DataFrame(0, index=consensus_peaks_df.index, columns=[sample])
    
    try:
        # Check if peakfile exists and is not empty
        if os.path.exists(peakfile_path) and os.path.getsize(peakfile_path) > 0:
            sample_peaks = bedtools.BedTool(peakfile_path)
            intersect_result = consensus_peaks.intersect(
                sample_peaks,
                g=chrom_file, 
                wa=True,
                c=True
            )
            
            if len(intersect_result) > 0:
                intersect_df = intersect_result.to_dataframe()
                # Check column structure - name should be in column 3 (0-indexed)
                if len(intersect_df.columns) >= 4:
                    intersect_df.columns = ['chrom', 'start', 'end', 'name', 'count'] if len(intersect_df.columns) == 5 else ['chrom', 'start', 'end', 'name']
                    if 'count' in intersect_df.columns:
                        intersect_df = intersect_df.set_index('name')[['count']]
                        intersect_df.columns = [sample]
                        result = intersect_df.reindex(consensus_peaks_df.index, fill_value=0)
                    else:
                        # If no count column, just mark as 1 if intersected
                        intersect_df = intersect_df.set_index('name')
                        result.loc[intersect_df.index, sample] = 1
        else:
            # Empty peakfile, all zeros
            result[sample] = 0
    except Exception as e:
        print(f"Error occurred while processing sample {sample}: {e}")
        import traceback
        traceback.print_exc()
    finally:
        result.T.to_csv(quant_support_path)