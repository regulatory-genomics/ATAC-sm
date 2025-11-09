
rule fastqc_1:
    input:
        lambda w: get_units_fastqs(w)[0],
    output:
        html=os.path.join(result_path,"report","{sample_run}_fastqc_1.html"),
        zip=os.path.join(result_path,"report","{sample_run}_fastqc_1.zip"), # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    params:
        extra="--quiet",
        mem_overhead_factor=0.1,
    log:
        "logs/rules/fastqc/{sample_run}.log",
    threads: 2
    resources:
        mem_mb = 1024,
    wrapper:
        "v7.6.1/bio/fastqc"

rule fastqc_2:
    input:
        lambda w: get_units_fastqs(w)[1],
    output:
        html=os.path.join(result_path,"report","{sample_run}_fastqc_2.html"),
        zip=os.path.join(result_path,"report","{sample_run}_fastqc_2.zip"), # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    params:
        extra="--quiet",
        mem_overhead_factor=0.1,
    log:
        "logs/rules/fastqc/{sample_run}.log",
    threads: 2
    resources:
        mem_mb = 1024,
    wrapper:
        "v7.6.1/bio/fastqc"

rule trim_galore_pe:
    input:
        get_units_fastqs,
    output:
        fasta_fwd=temp(os.path.join(result_path,"trimmed","{sample_run}_1.fq.gz")),
        report_fwd=temp(os.path.join(result_path,"trimmed","{sample_run}_1._trimming_report.txt")),
        fasta_rev=temp(os.path.join(result_path,"trimmed","{sample_run}_2.fq.gz")),
        report_rev=temp(os.path.join(result_path,"trimmed","{sample_run}_2._trimming_report.txt")),
    threads: 2
    log:
        "logs/trim_galore/{sample_run}.log"
    wrapper:
        "v7.6.0/bio/trim_galore/pe"


rule tss_coverage:
    input:
        bam = os.path.join(result_path,"results","{sample}","mapped","{sample}.filtered.bam"),
        bai = os.path.join(result_path,"results","{sample}","mapped","{sample}.filtered.bam.bai"),
        chromosome_sizes = config["chromosome_sizes"],
        unique_tss = config["unique_tss"],
    output:
        tss_hist = os.path.join(result_path,"results","{sample}","{sample}.tss_histogram.csv"),
    params:
        noise_upper = ( config["tss_slop"] * 2 ) - config["noise_lower"],
        double_slop = ( config["tss_slop"] * 2 ),
        genome_size = config["genome_size"],
        tss_slop = config["tss_slop"],
        noise_lower = config["noise_lower"],
    resources:
        mem_mb=config.get("mem", "16000"),
    threads: 4*config.get("threads", 2)
    conda:
        "../envs/pybedtools.yaml",
    log:
        "logs/rules/coverage_{sample}.log"
    shell:
        """
        echo "base,count" > {output.tss_hist};
        bedtools slop -b {params.tss_slop} -i {input.unique_tss} -g {input.chromosome_sizes} | \
            bedtools coverage -a - -b {input.bam} -d -sorted | \
            awk '{{if($6 == "+"){{ counts[$7] += $8;}} else counts[{params.double_slop} - $7 + 1] += $8;}} END {{ for(pos in counts) {{ if(pos < {params.noise_lower} || pos > {params.noise_upper}) {{ noise += counts[pos] }} }}; average_noise = noise /(2 * {params.noise_lower}); for(pos in counts) {{print pos-2000-1","(counts[pos]/average_noise) }} }}' | \
            sort -t "," -k1,1n >> {output.tss_hist} ;
        """


rule ataqv:
    input:
        bam = os.path.join(result_path, "results", "{sample}", "mapped", "{sample}.filtered.bam"),
        bai = os.path.join(result_path, "results", "{sample}", "mapped", "{sample}.filtered.bam.bai"),
        peak_file = os.path.join(result_path, "results", "{sample}", "peaks", "{sample}_peaks.narrowPeak"),
        tss_file = os.path.join(result_path, "genome", "tss.bed"),
        excl_regs_file = config["blacklisted_regions"],
        autosom_ref_file = os.path.join(result_path, "genome", "autosomes.txt")
    output:
        json = os.path.join(result_path, "results", "{sample}", "{sample}.ataqv.json")
    params:
        organism = config.get("organism", "hg38"),
        mito_name = config.get("mito_name", config.get("mitochondria_name", "chrM")),
        args = config.get("ataqv_args", "--ignore-read-groups"),
        prefix = "{sample}"
    log:
        "logs/ataqv/{sample}.log"
    conda:
        "/data2/litian/macrophage/script/atacseq_pipeline/.snakemake/conda/22ecf29b8590ba142cad5ef551bfa62e_"
    threads: config.get("ataqv_threads", 2)
    shell:
        """
        ataqv \
            {params.args} \
            --mitochondrial-reference-name {params.mito_name} \
            --peak-file {input.peak_file} \
            --tss-file {input.tss_file} \
            --excluded-region-file {input.excl_regs_file} \
            --autosomal-reference-file {input.autosom_ref_file} \
            --metrics-file {output.json} \
            --threads {threads} \
            --name {params.prefix} \
            {params.organism} \
            {input.bam} 2> {log}
        """

rule mkarv:
    input:
        jsons = expand(os.path.join(result_path, "results", "{sample}", "{sample}.ataqv.json"), sample=samples.keys())
    output:
        html = directory(os.path.join(result_path, "ataqv_report"))
    params:
        args = config.get("mkarv_args", "")
    log:
        "logs/mkarv/mkarv.log"
    conda:
        "/data2/litian/macrophage/script/atacseq_pipeline/.snakemake/conda/22ecf29b8590ba142cad5ef551bfa62e_"
    threads: config.get("mkarv_threads", 2)
    shell:
        """
        mkdir -p {output.html}
        mkarv \
            {params.args} \
            --concurrency {threads} \
            --force \
            {output.html}/ \
            {input.jsons} 2> {log}

        """

