
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
        
# generate consensus regions using (py)bedtools
rule get_consensus_regions:
    input:
        summits_bed = expand(os.path.join(result_path,"important_processed","peaks","{sample}_summits.bed"), sample=samples.keys()),
        blacklisted_regions = config["refs"]["blacklist"],
        chromosome_sizes = config["refs"]["chrom_sizes"],
    output:
        consensus_regions = os.path.join(result_path,"downstream_res","annotation","consensus_regions.bed"),
    params:
        slop_extension = config["peaks"]["slop_extension"],
    resources:
        mem_mb=config["resources"].get("mem_mb", 16000),
    threads: config["resources"].get("threads", 2)
    conda:
        "../envs/pybedtools.yaml",
    log:
        "logs/rules/annotation/get_consensus_regions.log"
    script:
        "../scripts/get_consensus_regions.py"

# quantify coverage based on consensus regions support for every sample
rule quantify_support_sample:
    input:
        consensus_regions = os.path.join(result_path,"downstream_res","annotation","consensus_regions.bed"),
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
        mem_mb=config["resources"].get("mem_mb", 32000),
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
        mem_mb=config["resources"].get("mem_mb", 16000),
    threads: config["resources"].get("threads", 2)
    conda:
        "../envs/pybedtools.yaml",
    log:
        "logs/rules/quantification/map_consensus_tss.log"
    script:
        "../scripts/map_consensus_tss.py"