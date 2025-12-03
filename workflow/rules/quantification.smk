
# create sample annotation file based on MultiQC general stats
rule sample_annotation:
    input:
        multiqc_stats = os.path.join(result_path, "report", "multiqc_report_data", "multiqc_general_stats.txt"),
    output:
        sample_annot = os.path.join(result_path, "downstream_res", "annotation", "sample_annotation.csv"),
    resources:
        mem_mb=config["resources"].get("mem_mb", 16000),
    threads: config["resources"].get("threads", 2)
    log:
        "logs/rules/annotation/sample_annotation.log"
    run:
        annot_df = pd.read_csv(input.multiqc_stats, delimiter='\t', index_col=0)
        # samples is a dict, convert to list of keys for indexing
        sample_list = list(samples.keys())
        annot_df = annot_df.loc[sample_list,:]
        annot_df.columns = [col.split("mqc-generalstats-")[1].replace("the_atac_seq_pipeline-", "").replace('-', '_') for col in annot_df.columns]
        annot_df.index.names = ['sample_name']
        annot_df.to_csv(output.sample_annot)

# generate promoter regions using (py)bedtools
rule get_promoter_regions:
    input:
        gencode_gtf = config["refs"]["gencode_gtf"],
        chromosome_sizes = config["refs"]["chrom_sizes"],
        genome_fasta = config["refs"]["fasta"],
    output:
        promoter_regions = os.path.join(result_path,"downstream_res","annotation","promoter_regions.bed"),
        promoter_annot = os.path.join(result_path,"downstream_res","annotation","promoter_annotation.csv"),
    params:
        proximal_size_up = config["annotation"]["promoter"]["up"],
        proximal_size_dn = config["annotation"]["promoter"]["down"],
    resources:
        mem_mb = config["resources"].get("mem_mb", 16000),
    threads: config["resources"].get("threads", 2)
    conda:
        "../envs/pybedtools.yaml",
    log:
        "logs/rules/annotation/get_promoter_regions.log"
    script:
        "../scripts/get_promoter_regions.py"

# quantify coverage based on consensus regions support for every sample
rule quantify_support_sample:
    input:
        consensus_regions = os.path.join(result_path,"downstream_res", "merged_peaks", "merged_peaks.bed"),
        peakfile = os.path.join(result_path,"important_processed","peaks", "{sample}_summits.bed"),
        chromosome_sizes = config["refs"]["chrom_sizes"],
    output:
        quant_support = os.path.join(result_path,"downstream_res","quantification","{sample}_quantification_support_counts.csv"),
    params:
        slop_extension = config["peaks"]["slop_extension"],
    resources:
        mem_mb=config["resources"].get("mem_mb", 16000),
    threads: config["resources"].get("threads", 2)
    conda:
        "../envs/pybedtools.yaml",
    log:
        "logs/rules/quantification/quantify_support_sample_{sample}.log"
    script:
        "../scripts/quantify_support_sample.py"

# quantify coverage based on consensus regions counts for every sample
rule quantify_counts_sample:
    input:
        regions = os.path.join(result_path,"downstream_res","annotation","{kind}_regions.bed"),
        bamfile = os.path.join(result_path,"important_processed","bam","{sample}.filtered.bam"),
        chromosome_sizes = config["refs"]["chrom_sizes"],
    output:
        quant_counts = os.path.join(result_path,"downstream_res","quantification","{sample}_quantification_{kind}_counts.csv"),
    resources:
        mem_mb=config["resources"].get("mem_mb", 16000),
    threads: config["resources"].get("threads", 2)
    conda:
        "../envs/pybedtools.yaml",
    log:
        "logs/rules/quantification/quantify_sample_{sample}_{kind}.log"
    script:
        "../scripts/quantify_counts_sample.py"
        
# aggregate quantification of counts/support of all samples
rule quantify_aggregate:
    input:
        get_quantifications,
    wildcard_constraints:
        kind="support|consensus|promoter|TSS"
    output:
        os.path.join(result_path,"downstream_res","quantification","{kind}_counts.csv"),
    resources:
        mem_mb=3*config["resources"].get("mem_mb", 32000),
    threads: 2*config["resources"].get("threads", 2)
    conda:
        "../envs/datamash.yaml",
    log:
        "logs/rules/quantification/quantify_aggregate_{kind}.log"
    shell:
        """
        awk 'NR==1 {{print; next}} FNR>1 {{print}}' {input} | datamash transpose -t ',' > {output}
        """
        
# aggregate HOMER motif enrichment results for all QC'd samples into one CSV
rule homer_aggregate:
    input:
        expand(os.path.join(result_path,"downstream_res","homer","{sample}","knownResults.txt"), sample=samples.keys()),
    output:
        os.path.join(result_path,"downstream_res","quantification","HOMER_knownMotifs.csv"),
    resources:
        mem_mb=config["resources"].get("mem_mb", 32000),
    threads: config["resources"].get("threads", 2)
    log:
        "logs/rules/quantification/homer_aggregate.log"
    run:
        combined_df = pd.DataFrame()

        for file_path in input:
            if os.path.getsize(file_path) > 0:
                sample = file_path.split('/')[-2]  # Extract sample name from path: .../homer/{sample}/knownResults.txt
                df = pd.read_csv(file_path, sep='\t')

                # replace white space in column names
                df.columns = [col.replace(' ', '_') for col in df.columns]

                # Remove columns that start with '#_' (unique per sample)
                df = df.loc[:, ~df.columns.str.startswith('#_')]

                # add sample name
                df.insert(0, 'sample_name', sample)

                if combined_df.shape[0]==0:
                    combined_df = df
                else:
                    combined_df = pd.concat([combined_df, df], ignore_index=True)

        combined_df.to_csv(output[0], index=False)

# map consensus regions to closest TSS per gene
rule map_consensus_tss:
    input:
        region_annotation = os.path.join(result_path,'downstream_res','annotation',"consensus_annotation.csv"),
        consensus_counts = os.path.join(result_path,"downstream_res","quantification","consensus_counts.csv"),
    output:
        tss_counts = os.path.join(result_path,"downstream_res","quantification","TSS_counts.csv"),
        tss_annot = os.path.join(result_path,"downstream_res","annotation","TSS_annotation.csv"),
        tss_bed = os.path.join(result_path,"downstream_res","annotation","TSS_regions.bed"),
    resources:
        # map_consensus_tss can be memory-heavy on large consensus matrices
        mem_mb=4*config["resources"].get("mem_mb", 16000),
        runtime = 240,
    threads: 4*config["resources"].get("threads", 2)
    conda:
        "../envs/pybedtools.yaml",
    log:
        "logs/rules/quantification/map_consensus_tss.log"
    script:
        "../scripts/map_consensus_tss.py"


rule merge_peaks:
    input:
        peak_calls = expand(
            os.path.join(result_path, "important_processed", "peaks", "{sample}_peaks.narrowPeak"),
            sample=get_samples_passing_qc(),
        )
    output:
        os.path.join(result_path, "downstream_res", "merged_peaks", "merged_peaks.bed"),
    threads:
        config["resources"].get("threads", 1)
    resources:
        mem_mb=40000,
        runtime = 30,
    params:
        chrom_sizes = config["refs"]["chrom_sizes"],
    conda:
        "../envs/mergepeak.yaml",
    log:
        "logs/rules/merge_peaks.log"
    shell:
        """
        merge_peaks --chrom-sizes {params.chrom_sizes} \
            --half-width 250 \
            --score-threshold 5 \
            --overlap-threshold 1 \
            --output {output} \
            {input.peak_calls}
        """

# count reads per peak for each sample using bedtools multicov
rule count_peaks_sample:
    input:
        merged_peaks = os.path.join(result_path, "downstream_res", "merged_peaks", "merged_peaks.bed"),
        bam = os.path.join(result_path, "important_processed", "bam", "{sample}.filtered.bam"),
        bai = os.path.join(result_path, "important_processed", "bam", "{sample}.filtered.bam.bai"),
    output:
        counts = os.path.join(result_path, "downstream_res", "quantification", "merged_peaks", "{sample}_counts.txt"),
    resources:
        mem_mb=3*config["resources"].get("mem_mb", 16000),
        runtime = 200,
    threads: 4*config["resources"].get("threads", 2)
    conda:
        "../envs/pybedtools.yaml",
    log:
        "logs/rules/quantification/count_peaks_sample_{sample}.log"
    wildcard_constraints:
        sample="|".join(samples.keys())
    shell:
        """
        mkdir -p $(dirname {output.counts})
        # Run bedtools multicov and extract only the count column (last column)
        bedtools multicov -bams {input.bam} -bed {input.merged_peaks} | awk '{{print $NF}}' > {output.counts}
        # Add sample name as header
        sed -i "1s/^/{wildcards.sample}\\n/" {output.counts}
        """

# aggregate count matrix from all samples
rule count_peaks_matrix:
    input:
        merged_peaks = os.path.join(result_path, "downstream_res", "merged_peaks", "merged_peaks.bed"),
        sample_counts = expand(
            os.path.join(result_path, "downstream_res", "quantification", "merged_peaks", "{sample}_counts.txt"),
            sample=samples.keys()
        ),
        sample_annotation = annotation_sheet_path,
    output:
        count_matrix = os.path.join(result_path, "downstream_res", "quantification", "merged_peaks_count_matrix.txt"),
    resources:
        mem_mb=config["resources"].get("mem_mb", 32000),
    threads: config["resources"].get("threads", 2)
    log:
        "logs/rules/quantification/count_peaks_matrix.log"
    run:
        import pandas as pd
        import os
        
        # Read merged peaks file to determine number of columns
        peaks_df = pd.read_csv(input.merged_peaks, sep='\t', header=None, nrows=1)
        n_cols = peaks_df.shape[1]
        
        # Read full peaks file with appropriate column names
        if n_cols >= 4:
            peaks_df = pd.read_csv(input.merged_peaks, sep='\t', header=None, 
                                  names=['chr', 'start', 'end', 'peak_id'])
        else:
            peaks_df = pd.read_csv(input.merged_peaks, sep='\t', header=None, 
                                  names=['chr', 'start', 'end'])
            # Create peak_id column if it doesn't exist
            peaks_df['peak_id'] = peaks_df.index.to_series().apply(lambda x: f"peak_{x+1:06d}")
        
        # Extract coordinate columns
        final_matrix = peaks_df[['chr', 'start', 'end', 'peak_id']].copy()
        
        # Read annotation file to maintain sample order (as in original bash script)
        annot_df = pd.read_csv(input.sample_annotation)
        # Get unique sample names in order (maintaining first occurrence order)
        sample_order = annot_df['sample_name'].drop_duplicates().tolist()
        # Filter to only samples that exist in samples.keys()
        sample_order = [s for s in sample_order if s in samples.keys()]
        
        # Read count files in the correct order and add as columns
        n_peaks = final_matrix.shape[0]
        for sample in sample_order:
            count_file = os.path.join(result_path, "downstream_res", "quantification", "merged_peaks", f"{sample}_counts.txt")
            if os.path.exists(count_file):
                # Skip first line (header with sample name) and read counts as a single column
                counts = pd.read_csv(count_file, header=None, skiprows=1, names=[sample])
                # Basic sanity check: number of rows must match number of peaks
                if counts.shape[0] != n_peaks:
                    raise ValueError(f"Count file {count_file} has {counts.shape[0]} rows, "
                                     f"but merged_peaks has {n_peaks} rows.")
                # Assign counts as a new column; avoid concat to dodge MultiIndex/axis-union issues
                final_matrix[sample] = counts[sample].to_numpy()
        
        # Save final matrix
        os.makedirs(os.path.dirname(output.count_matrix), exist_ok=True)
        final_matrix.to_csv(output.count_matrix, sep='\t', index=False)
        
        # Cleanup temporary count files
        for count_file in input.sample_counts:
            if os.path.exists(count_file):
                os.remove(count_file)