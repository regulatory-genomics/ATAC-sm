#!/bin/env python

#### libraries
import pandas as pd
import pybedtools as bedtools

# map region to gene and classify if TSS
def map_region(x):
    tmp_dist = x.loc[x['homer_Distance_to_TSS'].abs().idxmin(),'homer_Distance_to_TSS']
    if (tmp_dist>=TSS_up) and (tmp_dist<=TSS_dn):
        return x.loc[x['homer_Distance_to_TSS'].abs().idxmin(),:]
    else:
        return None

#### configurations

# input
region_annotation_path = snakemake.input["region_annotation"]
consensus_counts_path = snakemake.input["consensus_counts"]

# output
tss_counts_path = snakemake.output["tss_counts"]
tss_annot_path = snakemake.output["tss_annot"]
tss_bed_path = snakemake.output["tss_bed"]

# parameters
TSS_up = -snakemake.config["annotation"]["promoter"]["up"]
TSS_dn = snakemake.config["annotation"]["promoter"]["down"]

# load annotations and consensus counts
try:
    annot_regions = pd.read_csv(region_annotation_path, index_col='peak_id')
except (KeyError, ValueError):
    # If peak_id is not the index, read without index and set it
    annot_regions = pd.read_csv(region_annotation_path)
    if 'peak_id' in annot_regions.columns:
        annot_regions.set_index('peak_id', inplace=True)
    else:
        # No peak_id column, create empty outputs
        pd.DataFrame().to_csv(tss_counts_path)
        pd.DataFrame().to_csv(tss_annot_path)
        with open(tss_bed_path, 'w') as f:
            pass
        exit(0)

consensus_counts = pd.read_csv(consensus_counts_path, index_col=0)

# Check if DataFrames are empty
if len(annot_regions) == 0 or len(consensus_counts) == 0:
    # Create empty output files
    pd.DataFrame().to_csv(tss_counts_path)
    pd.DataFrame().to_csv(tss_annot_path)
    # Create empty BED file
    with open(tss_bed_path, 'w') as f:
        pass
else:
    # Check if required columns exist
    if 'homer_Nearest_Ensembl' not in annot_regions.columns or 'homer_Distance_to_TSS' not in annot_regions.columns:
        # Missing required columns, create empty outputs
        pd.DataFrame().to_csv(tss_counts_path)
        pd.DataFrame().to_csv(tss_annot_path)
        with open(tss_bed_path, 'w') as f:
            pass
    else:
        # Reset index to make peak_id a column for groupby operations
        annot_regions_reset = annot_regions.reset_index()
        
        # Filter out rows where homer_Nearest_Ensembl is NaN
        annot_regions_reset = annot_regions_reset[annot_regions_reset['homer_Nearest_Ensembl'].notna()]
        
        if len(annot_regions_reset) == 0:
            # No valid Ensembl IDs, create empty outputs
            pd.DataFrame().to_csv(tss_counts_path)
            pd.DataFrame().to_csv(tss_annot_path)
            with open(tss_bed_path, 'w') as f:
                pass
        else:
            # map regions to closest regions and subset for regions in TSS proximity
            TSS_regions = annot_regions_reset.groupby('homer_Nearest_Ensembl').apply(map_region)
            TSS_regions = TSS_regions.dropna(axis=0, how='all')
            
            if len(TSS_regions) == 0:
                # No TSS regions found, create empty outputs
                pd.DataFrame().to_csv(tss_counts_path)
                pd.DataFrame().to_csv(tss_annot_path)
                with open(tss_bed_path, 'w') as f:
                    pass
            else:
                # After groupby().apply(), the result has the groupby key (gene name) as index
                # and peak_id should be in the columns (from the returned Series/DataFrame)
                # Reset index to make gene names a column and peak_id accessible
                TSS_regions = TSS_regions.reset_index()
                
                # Ensure peak_id column exists
                if 'peak_id' not in TSS_regions.columns:
                    # peak_id should be in the data returned by map_region
                    # If not found, we can't proceed
                    pd.DataFrame().to_csv(tss_counts_path)
                    pd.DataFrame().to_csv(tss_annot_path)
                    with open(tss_bed_path, 'w') as f:
                        pass
                else:
                    # Get peak_ids and gene names (from the groupby key, now in 'homer_Nearest_Ensembl' column after reset_index)
                    peak_ids = TSS_regions['peak_id'].values
                    # The gene names are in the index column after reset_index (the groupby key)
                    # The column name should be 'homer_Nearest_Ensembl' (the groupby column)
                    gene_names = TSS_regions['homer_Nearest_Ensembl'].values if 'homer_Nearest_Ensembl' in TSS_regions.columns else TSS_regions.index.values
                    
                    # Create mapping from peak_id to gene name
                    peak_to_gene = dict(zip(peak_ids, gene_names))
                    
                    # Only select peak_ids that exist in consensus_counts
                    valid_peak_ids = [pid for pid in peak_ids if pid in consensus_counts.index]
                    if len(valid_peak_ids) == 0:
                        TSS_counts = pd.DataFrame()
                    else:
                        TSS_counts = consensus_counts.loc[valid_peak_ids, :]
                        # Set index to gene names
                        TSS_counts.index = [peak_to_gene.get(pid, pid) for pid in TSS_counts.index]
                    TSS_counts.to_csv(tss_counts_path)
                    
                    # subset the consensus annotation by the successfully mapped consenesus regions
                    if len(valid_peak_ids) == 0:
                        TSS_annot = pd.DataFrame()
                    else:
                        TSS_annot = annot_regions.loc[valid_peak_ids, :].copy()
                        TSS_annot.reset_index(inplace=True)
                        # Map peak_ids to gene names
                        TSS_annot.index = [peak_to_gene.get(pid, pid) for pid in TSS_annot['peak_id']]
                    TSS_annot.to_csv(tss_annot_path)
                    
                    # save bed file of TSS regions
                    if len(TSS_annot) == 0 or 'gencode_chr' not in TSS_annot.columns:
                        # Create empty BED file
                        with open(tss_bed_path, 'w') as f:
                            pass
                    else:
                        TSS_bed_df = TSS_annot.sort_values(by="peak_id")[["gencode_chr", 'gencode_start', 'gencode_end', 'homer_Nearest_Ensembl']]
                        TSS_bed = bedtools.BedTool.from_dataframe(TSS_bed_df)
                        TSS_bed.saveas(tss_bed_path)