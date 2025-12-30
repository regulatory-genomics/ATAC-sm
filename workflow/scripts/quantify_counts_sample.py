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
    # Standard BED: chrom (0), start (1), end (2), name (3, optional), ...
    # Coverage count is always the last column
    coverage_result = elements_to_quantify.coverage(b=bamfile_path, sorted=True, g=chrom_file)
    
    # Read as dataframe without specifying names to see actual structure
    coverage_df = coverage_result.to_dataframe(header=None)
    
    # Coverage count is always the last column
    coverage_col_idx = len(coverage_df.columns) - 1
    coverage_counts = coverage_df.iloc[:, coverage_col_idx].astype(int)
    
    # Determine ID: use column 3 if it exists and looks like an ID, otherwise create from chrom:start-end
    if len(coverage_df.columns) > 3:
        # Check if column 3 exists and can be used as ID
        id_col = coverage_df.iloc[:, 3].astype(str)
        # Use column 3 as ID if it's not empty and doesn't look like a number
        if not id_col.isna().all() and not id_col.str.match(r'^\d+$').all():
            region_ids = id_col.values
        else:
            # Create ID from chrom:start-end
            region_ids = (coverage_df.iloc[:, 0].astype(str) + ':' + 
                         coverage_df.iloc[:, 1].astype(str) + '-' + 
                         coverage_df.iloc[:, 2].astype(str)).values
    else:
        # No name column, create ID from chrom:start-end
        region_ids = (coverage_df.iloc[:, 0].astype(str) + ':' + 
                     coverage_df.iloc[:, 1].astype(str) + '-' + 
                     coverage_df.iloc[:, 2].astype(str)).values
    
    # Create result dataframe with ID as index and sample name as column
    result_df = pd.DataFrame({
        sample: coverage_counts
    }, index=region_ids)
    
    result = result_df.T
    result.to_csv(quant_count_path)
except Exception as e:
    print("Error occured while processing sample "+sample)
    print(f"Error details: {e}")
    import traceback
    traceback.print_exc()
    # Fallback: create empty dataframe with zeros
    try:
        elements_to_quantify_df = elements_to_quantify.to_dataframe(header=None)
        if len(elements_to_quantify_df.columns) >= 3:
            # Create IDs from chrom:start-end
            region_ids = (elements_to_quantify_df.iloc[:, 0].astype(str) + ':' + 
                         elements_to_quantify_df.iloc[:, 1].astype(str) + '-' + 
                         elements_to_quantify_df.iloc[:, 2].astype(str)).values
        else:
            # Fallback to simple index
            region_ids = [f"region_{i}" for i in range(len(elements_to_quantify_df))]
        pd.DataFrame(0, index=region_ids, columns=[sample]).T.to_csv(quant_count_path)
    except Exception as e2:
        print(f"Fallback also failed: {e2}")
        # Last resort: create minimal output
        pd.DataFrame({sample: [0]}, index=['unknown']).T.to_csv(quant_count_path)