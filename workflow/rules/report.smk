
aligner = config["alignment"].get("tool", "bowtie2")
_ALIGNER_LOG_SUFFIX = "txt" if aligner == "bowtie2" else "bwa.log"


rule collect_align_stats:
    input:
        expand(os.path.join(result_path, 'results', '{sample_run}', '{sample_run}.align.stats.tsv'), 
               sample_run=annot.index.tolist())
    output:
        report_tsv = os.path.join(result_path, 'report', 'align_stats_report.tsv')
    resources:
        mem_mb=config["resources"].get("mem_mb", 1000),
    threads: config["resources"].get("threads", 1)
    log:
        os.path.join("logs", "rules", "collect_align_stats.log")
    shell:
        """
        echo -e "sample_run\\tmetric\\tvalue" > {output.report_tsv}
        for stats_file in {input}; do
            if [ -f "$stats_file" ] && [ -s "$stats_file" ]; then
                sample_run=$(basename "$stats_file" .align.stats.tsv)
                while IFS=$'\\t' read -r metric value || [ -n "$metric" ]; do
                    # Skip empty lines
                    [ -z "$metric" ] && continue
                    echo -e "$sample_run\\t$metric\\t$value" >> {output.report_tsv}
                done < "$stats_file"
            fi
        done
        """


rule symlink_run_stats:
    input:
        align_stats = os.path.join(result_path, 'results', "{sample_run}", '{sample_run}.align.stats.tsv'),
        mapped_log = lambda w: os.path.join(result_path, 'results', w.sample_run, 'mapped', f'{w.sample_run}.{_ALIGNER_LOG_SUFFIX}'),
        samblaster_log = os.path.join(result_path, 'results', "{sample_run}", 'mapped', '{sample_run}.samblaster.log'),
        flagstat_log = os.path.join(result_path, 'results', "{sample_run}", 'mapped', '{sample_run}.samtools_flagstat.log'),
    wildcard_constraints:
        sample_run="|".join(annot.index.tolist())
    output:
        align_stats = os.path.join(result_path, 'report', '{sample_run}.align.stats.tsv'),
        mapped_log = os.path.join(result_path, 'report', '{sample_run}.' + _ALIGNER_LOG_SUFFIX),
        samblaster_log = os.path.join(result_path, 'report', '{sample_run}.samblaster.log'),
        flagstat_log = os.path.join(result_path, 'report', '{sample_run}.samtools_flagstat.log'),
    resources:
        mem_mb=config["resources"].get("mem_mb", 1000),
    threads: config["resources"].get("threads", 1)
    log:
        os.path.join("logs", "rules", "symlink_run_stats_{sample_run}.log")
    shell:
        """
        ln -sfn $(realpath --relative-to=$(dirname {output.align_stats}) {input.align_stats}) {output.align_stats}
        ln -sfn $(realpath --relative-to=$(dirname {output.mapped_log}) {input.mapped_log}) {output.mapped_log}
        ln -sfn $(realpath --relative-to=$(dirname {output.samblaster_log}) {input.samblaster_log}) {output.samblaster_log}
        ln -sfn $(realpath --relative-to=$(dirname {output.flagstat_log}) {input.flagstat_log}) {output.flagstat_log}
        """

rule symlink_stats:
    input:
        stats_tsv = os.path.join(result_path, 'results', "{sample}", '{sample}.stats.tsv'),
        tss_csv = os.path.join(result_path, 'results', "{sample}", '{sample}.tss_histogram.csv'),
        macs2_log = os.path.join(result_path, 'results', "{sample}", 'peaks', '{sample}.macs2.log'),
        peaks_xls = os.path.join(result_path, 'results', "{sample}", 'peaks', '{sample}_peaks.xls'),
    output:
        stats_tsv = os.path.join(result_path, 'report', '{sample}.stats.tsv'),
        tss_csv = os.path.join(result_path, 'report', '{sample}_TSS.csv'),
        macs2_log = os.path.join(result_path, 'report', '{sample}.macs2.log'),
        peaks_xls = os.path.join(result_path, 'report', '{sample}_peaks.xls'),
    resources:
        mem_mb=config["resources"].get("mem_mb", 1000),
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

rule multiqc:
    input:
        expand(os.path.join(result_path,"bam","{sample}", "{sample}.filtered.bam"), sample=samples.keys()),
        expand(os.path.join(result_path,"results","{sample}","peaks","{sample}_peaks.narrowPeak"), sample=samples.keys()),
        expand(os.path.join(result_path, 'report', '{sample}_peaks.xls'), sample=samples.keys()), # representing symlinked stats
        # Collect fastqc files from all runs (sample_run format)
        expand(os.path.join(result_path, 'report', '{sample_run}_fastqc_1.html'), sample_run=annot.index.tolist()),
        expand(os.path.join(result_path, 'report', '{sample_run}_fastqc_2.html'), sample_run=annot.index.tolist()),
        # Collect per-run stats and logs for MultiQC
        expand(os.path.join(result_path, 'report', '{sample_run}.align.stats.tsv'), sample_run=annot.index.tolist()),
        expand(os.path.join(result_path, 'report', '{sample_run}.' + _ALIGNER_LOG_SUFFIX), sample_run=annot.index.tolist()),
        expand(os.path.join(result_path, 'report', '{sample_run}.samblaster.log'), sample_run=annot.index.tolist()),
        expand(os.path.join(result_path, 'report', '{sample_run}.samtools_flagstat.log'), sample_run=annot.index.tolist()),
        sample_annotation = config["project"]["samples"],
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
        multiqc_configs = "{{'title': '{name}', 'intro_text': 'Quality Control Metrics of the ATAC-seq pipeline.', 'annotation': '{annot}', 'genome': '{genome}', 'exploratory_columns': {exploratory_columns}, 'skip_versions_section': true,'custom_content': {custom_content}}}".format(name = config["project"]["name"], annot = config["project"]["samples"], genome = config["project"]["genome"], exploratory_columns = config["project"].get("annot_columns", "[]"), custom_content = config.get("custom_content", "")),
    resources:
        mem_mb=config["resources"].get("mem_mb", 16000),
    threads: config["resources"].get("threads", 2)
    conda:
        "../envs/multiqc.yaml",
    log:
        "logs/rules/multiqc.log"
    shell:
        """
        multiqc {params.result_path}/report --force --verbose --outdir {params.result_path}/report --filename multiqc_report.html --cl-config "{params.multiqc_configs}"
        """

# visualize sample annotation (including QC metrics)
rule plot_sample_annotation:
    input:
        sample_annotation = config["project"]["samples"],
        sample_annotation_w_QC = os.path.join(result_path, "counts", "sample_annotation.csv"),
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