#!/bin/env python

#### libraries
import pandas as pd
import pybedtools as bedtools

#### configurations

# input
regions_path = snakemake.input["regions"]
bamfile_path = snakemake.input["bamfile"]
chrom_file = snakemake.input["chromosome_sizes"]

# output
quant_count_path = snakemake.output["quant_counts"]

# parameters
sample = snakemake.wildcards["sample"]

# Load and sort BED file to match genome file order (required for sorted=True)
elements_to_quantify = bedtools.BedTool(regions_path).sort(g=chrom_file)

print("Processing "+sample)
try:
    # bedtools coverage output: preserves all BED columns + adds coverage count as last column
    # ID is typically in column 3 (0-indexed), coverage count is the last column
    coverage_result = elements_to_quantify.coverage(b=bamfile_path, sorted=True, g=chrom_file)
    
    # Read as dataframe without specifying names to see actual structure
    coverage_df = coverage_result.to_dataframe(header=None)
    
    # ID column is column 3 (index 3), coverage count is the last column
    id_col_idx = 3
    coverage_col_idx = len(coverage_df.columns) - 1
    
    # Create result dataframe with ID as index and sample name as column
    result_df = pd.DataFrame({
        sample: coverage_df.iloc[:, coverage_col_idx].astype(int)
    }, index=coverage_df.iloc[:, id_col_idx])
    
    result = result_df.T
    result.to_csv(quant_count_path)
except Exception as e:
    print("Error occured while processing sample "+sample)
    print(f"Error details: {e}")
    # Fallback: create empty dataframe with zeros
    elements_to_quantify_df = elements_to_quantify.to_dataframe(
        names=["CHR", "START", "END", "ID"],
        usecols=['ID'],
        index_col='ID'
    )
    pd.DataFrame(0, index=elements_to_quantify_df.index, columns=[sample]).T.to_csv(quant_count_path)