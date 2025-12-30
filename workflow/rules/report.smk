
aligner = config["alignment"].get("tool", "bowtie2")
_ALIGNER_LOG_SUFFIX = "txt" if aligner == "bowtie2" else "bwa.log"

# Check if prealignments are enabled
prealign_enabled = config["alignment"].get("prealign", {}).get("enabled", True)
prealignments = config["alignment"].get("prealign", {}).get("indices", []) or []
has_prealignments = prealign_enabled and len(prealignments) > 0


rule collect_align_stats:
    input:
        expand(os.path.join(result_path, "report", "align_stats", '{sample}.align.stats.tsv'), 
                sample=samples.keys())
    output:
        report_tsv = os.path.join(result_path, 'report', 'align_stats_report.tsv')
    resources:
        mem_mb=config["resources"].get("mem_mb", 1000),
        runtime = 10
    threads: config["resources"].get("threads", 1)
    log:
        os.path.join("logs", "rules", "collect_align_stats.log")
    shell:
        """
        echo -e "sample\\tmetric\\tvalue" > {output.report_tsv}
        for stats_file in {input}; do
            if [ -f "$stats_file" ] && [ -s "$stats_file" ]; then
                sample=$(basename "$stats_file" .align.stats.tsv)
                while IFS=$'\\t' read -r metric value || [ -n "$metric" ]; do
                    # Skip empty lines
                    [ -z "$metric" ] && continue
                    echo -e "$sample\\t$metric\\t$value" >> {output.report_tsv}
                done < "$stats_file"
            fi
        done
        """


rule symlink_sample_stats:
    input:
        align_stats = os.path.join(result_path, "report", "align_stats", "{sample}.align.stats.tsv"),
        mapped_log = lambda w: os.path.join(result_path, 'logs',"align", w.sample, f'{w.sample}.{_ALIGNER_LOG_SUFFIX}') if aligner == "bwa-mem2" else os.path.join(result_path, 'logs',"align", f'{w.sample}.{_ALIGNER_LOG_SUFFIX}'),
        samblaster_log = lambda w: os.path.join(result_path, 'logs', 'align', w.sample, f'{w.sample}.samblaster.log') if aligner == "bwa-mem2" else os.path.join(result_path, 'logs', 'align', f'{w.sample}.samblaster.log'),
        flagstat_log = os.path.join(result_path, 'logs', 'align', '{sample}.samtools_flagstat.log'),
    wildcard_constraints:
        sample="|".join(samples.keys())
    output:
        align_stats = os.path.join(result_path, 'report', '{sample}.align.stats.tsv'),
        mapped_log = os.path.join(result_path, 'report', '{sample}.' + _ALIGNER_LOG_SUFFIX),
        samblaster_log = os.path.join(result_path, 'report', '{sample}.samblaster.log'),
        flagstat_log = os.path.join(result_path, 'report', '{sample}.samtools_flagstat.log'),
    resources:
        mem_mb=config["resources"].get("mem_mb", 1000),
        runtime = 10,
    threads: config["resources"].get("threads", 1)
    log:
        os.path.join("logs", "rules", "symlink_sample_stats_{sample}.log")
    shell:
        """
        ln -sfn $(realpath --relative-to=$(dirname {output.align_stats}) {input.align_stats}) {output.align_stats}
        ln -sfn $(realpath --relative-to=$(dirname {output.mapped_log}) {input.mapped_log}) {output.mapped_log}
        ln -sfn $(realpath --relative-to=$(dirname {output.samblaster_log}) {input.samblaster_log}) {output.samblaster_log}
        ln -sfn $(realpath --relative-to=$(dirname {output.flagstat_log}) {input.flagstat_log}) {output.flagstat_log}
        """

rule symlink_stats:
    input:
        stats_tsv = os.path.join(result_path, 'report', "peaks_stats", '{sample}.stats.tsv'),
        tss_csv = os.path.join(result_path, 'report', "tss_coverage", '{sample}.tss_histogram.csv'),
        macs2_log = os.path.join(result_path, 'important_processed', 'peaks', '{sample}.macs2.log'),
        peaks_xls = os.path.join(result_path, 'important_processed', 'peaks', '{sample}_peaks.xls'),
    wildcard_constraints:
        sample="|".join(samples.keys())
    output:
        stats_tsv = os.path.join(result_path, 'report', '{sample}.stats.tsv'),
        tss_csv = os.path.join(result_path, 'report', '{sample}_TSS.csv'),
        macs2_log = os.path.join(result_path, 'report', '{sample}.macs2.log'),
        peaks_xls = os.path.join(result_path, 'report', '{sample}_peaks.xls'),
    resources:
        mem_mb=config["resources"].get("mem_mb", 1000),
        runtime = 10,
    threads: config["resources"].get("threads", 1)
    log:
        os.path.join("logs", "rules", "symlink_stats_{sample}.log")
    shell:
        """
        ln -sfn $(realpath --relative-to=$(dirname {output.stats_tsv}) {input.stats_tsv}) {output.stats_tsv}
        ln -sfn $(realpath --relative-to=$(dirname {output.tss_csv}) {input.tss_csv}) {output.tss_csv}
        ln -sfn $(realpath --relative-to=$(dirname {output.macs2_log}) {input.macs2_log}) {output.macs2_log}
        ln -sfn $(realpath --relative-to=$(dirname {output.peaks_xls}) {input.peaks_xls}) {output.peaks_xls}
        """

rule generate_multiqc_sample_names:
    input:
        sample_annotation = annotation_sheet_path,
    output:
        sample_names_file = os.path.join(result_path, 'report', 'multiqc_sample_names.txt'),
    resources:
        mem_mb=config["resources"].get("mem_mb", 1000),
        runtime = 5,
    threads: 1
    log:
        "logs/rules/generate_multiqc_sample_names.log"
    run:
        import pandas as pd
        import os
        
        # Read annotation file
        annot_df = pd.read_csv(input.sample_annotation)
        
        # Create mapping file: fastq_filename_base -> sample_name
        # MultiQC extracts sample names from fastq filenames in the fastp JSON
        # We map those to the sample_name from annotation CSV
        # Use a dict to avoid duplicate mappings
        name_mapping = {}
        
        for _, row in annot_df.iterrows():
            sample_name = str(row['sample_name'])
            
            # Extract R1 filename base (what MultiQC will use)
            if 'R1' in row and pd.notna(row['R1']):
                r1_path = str(row['R1'])
                # Handle PEP derive modifier (raw_data|filename)
                if '|' in r1_path:
                    r1_path = r1_path.split('|', 1)[1]
                # Extract base filename without extension
                r1_basename = os.path.basename(r1_path)
                # Remove extensions (.fastq.gz, .fq.gz, etc.)
                r1_base = r1_basename.split('.')[0]
                # Map R1 base name -> sample_name
                name_mapping[r1_base] = sample_name
                
                # Also map R2 if it exists (for completeness)
                if 'R2' in row and pd.notna(row['R2']):
                    r2_path = str(row['R2'])
                    if '|' in r2_path:
                        r2_path = r2_path.split('|', 1)[1]
                    r2_basename = os.path.basename(r2_path)
                    r2_base = r2_basename.split('.')[0]
                    name_mapping[r2_base] = sample_name
        
        # Write the mapping file (tab-separated: original_name -> new_name)
        with open(output.sample_names_file, 'w') as f:
            for original_name, new_name in sorted(name_mapping.items()):
                f.write(f"{original_name}\t{new_name}\n")


# Rule to convert reproducibility QC JSON files to MultiQC custom content
rule reproducibility_to_multiqc:
    input:
        json_files = expand(
            os.path.join(result_path, "report", "reproducibility", "{group}.reproducibility.qc.json"),
            group=replicate_samples
        ) if len(replicate_samples) > 0 else [],
    output:
        multiqc_yaml = os.path.join(result_path, "report", "reproducibility", "reproducibility_mqc.yaml"),
        stats_tsv = os.path.join(result_path, "report", "reproducibility", "reproducibility_stats_mqc.tsv"),
    conda:
        "../envs/reproducibility.yaml"
    log:
        "logs/rules/reproducibility_to_multiqc.log"
    shell:
        """
        # If there are any reproducibility JSONs, run the conversion script
        if [ -n "{input.json_files}" ]; then
            python workflow/scripts/reproducibility_to_multiqc.py \
                {input.json_files} \
                --output-yaml {output.multiqc_yaml} \
                --output-stats {output.stats_tsv} \
                2>> {log}
        else
            # Create empty files if no replicate groups
            mkdir -p $(dirname {output.multiqc_yaml})
            touch {output.multiqc_yaml}
            touch {output.stats_tsv}
            echo "No replicate groups found, skipping reproducibility MultiQC conversion" >> {log}
        fi
        """

rule multiqc:
    input:
        expand(os.path.join(result_path,"important_processed","bam","{sample}.filtered.bam"), sample=samples.keys()),
        expand(os.path.join(result_path,"important_processed","peaks","{sample}_peaks.narrowPeak"), sample=samples.keys()),
        expand(os.path.join(result_path, 'report', '{sample}_peaks.xls'), sample=samples.keys()), # representing symlinked stats
        # collect fastp report from all runs (sample_run format)
        expand(os.path.join(result_path, 'report', 'fastp', '{sample_run}_fastp.html'), sample_run=annot.index.tolist()),
        expand(os.path.join(result_path, 'report', 'fastp', '{sample_run}_fastp.json'), sample_run=annot.index.tolist()),
        # Collect per-sample stats and logs for MultiQC
        expand(os.path.join(result_path, 'report', '{sample}.align.stats.tsv'), sample=samples.keys()),
        expand(os.path.join(result_path, 'report', '{sample}.' + _ALIGNER_LOG_SUFFIX), sample=samples.keys()),
        expand(os.path.join(result_path, 'report', '{sample}.samblaster.log'), sample=samples.keys()),
        expand(os.path.join(result_path, 'report', '{sample}.samtools_flagstat.log'), sample=samples.keys()),
        sample_annotation = annotation_sheet_path,
        sample_names_file = os.path.join(result_path, 'report', 'multiqc_sample_names.txt'),
        # Collect prealign stats if enabled (per sample_run, will be averaged per sample in MultiQC)
        prealign_stats = expand(os.path.join(result_path, 'report',"prealigned", '{sample_run}.prealign.stats.tsv'), sample_run=annot.index.tolist()) if has_prealignments else [],
        # Reproducibility QC for MultiQC (if replicate groups exist)
        reproducibility_yaml = os.path.join(result_path, "report", "reproducibility", "reproducibility_mqc.yaml") if len(replicate_samples) > 0 else [],
        reproducibility_stats = os.path.join(result_path, "report", "reproducibility", "reproducibility_stats_mqc.tsv") if len(replicate_samples) > 0 else [],
    output:
        multiqc_report = report(os.path.join(result_path,"report","multiqc_report.html"),
                                caption="../report/multiqc.rst",
                                category="{}_{}".format(config["project"]["name"], module_name),
                                subcategory="QC",
                                labels={
                                    "name": "MultiQC report",
                                    "type": "HTML",
                                }),
        multiqc_stats = os.path.join(result_path, "report", "multiqc_report_data", "multiqc_general_stats.txt"),
    params:
        result_path = result_path,
        multiqc_configs = "{{'title': '{name}', 'intro_text': 'Quality Control Metrics of the ATAC-seq pipeline.', 'annotation': '{annot}', 'genome': '{genome}', 'exploratory_columns': {exploratory_columns}, 'skip_versions_section': true}}".format(name = config["project"]["name"], annot = annotation_sheet_path, genome = config["project"]["genome"], exploratory_columns = config["project"].get("annot_columns", "[]"))
    resources:
        mem_mb=config["resources"].get("mem_mb", 16000),
        runtime = 20,
    threads: config["resources"].get("threads", 2)
    conda:
        "../envs/multiqc.yaml",
    log:
        "logs/rules/multiqc.log"
    shell:
        """
        multiqc {params.result_path}/report --force --verbose --outdir {params.result_path}/report --filename multiqc_report.html --replace-names {input.sample_names_file} --cl-config "{params.multiqc_configs}"
        """

# visualize sample annotation (including QC metrics)
rule plot_sample_annotation:
    input:
        sample_annotation = annotation_sheet_path,
        sample_annotation_w_QC = os.path.join(result_path, "downstream_res", "annotation", "sample_annotation.csv"),
    output:
        sample_annotation_plot = os.path.join(result_path,"report","sample_annotation.png"),
        sample_annotation_html = report(os.path.join(result_path,"report","sample_annotation.html"),
                       caption="../report/sample_annotation.rst",
                       category="{}_{}".format(config["project"]["name"], module_name),
                       subcategory="QC",
                       labels={
                           "name": "Sample annotation",
                           "type": "HTML",
                           }),
    log:
        "logs/rules/plot_sample_annotation.log",
    resources:
        mem_mb="4000",
    conda:
        "../envs/ggplot.yaml"
    script:
        "../scripts/plot_sample_annotation.R"
