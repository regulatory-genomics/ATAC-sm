
rule fastp:
    input:
        r1 = lambda w: get_units_fastqs(w)[0],
        r2 = lambda w: get_units_fastqs(w)[1],
    output:
        r1 = temp(os.path.join(result_path,"middle_files","trimmed","{sample_run}_1.fq.gz")),
        r2 = temp(os.path.join(result_path,"middle_files","trimmed","{sample_run}_2.fq.gz")),
        report_html = os.path.join(result_path,"report","fastp","{sample_run}_fastp.html"),
        report_json = os.path.join(result_path,"report","fastp","{sample_run}_fastp.json"),
    conda: 
        "../envs/fastp.yaml"
    resources:
        mem_mb=16000,
        runtime = 60,
    log:
        "logs/rules/fastp/{sample_run}.fastp.json"
    threads: 4
    shell:
        "fastp -i {input.r1} -I {input.r2} -o {output.r1} -O {output.r2} --detect_adapter_for_pe --trim_poly_g  --thread {threads} -j {output.report_json} -h {output.report_html}"


rule tss_coverage:
    input:
        bam = os.path.join(result_path,"important_processed","bam","{sample}.filtered.bam"),
        bai = os.path.join(result_path,"important_processed","bam","{sample}.filtered.bam.bai"),
        chromosome_sizes = config["refs"]["chrom_sizes"],
        unique_tss = config["refs"]["unique_tss"],
    output:
        tss_hist = os.path.join(result_path,"report","tss_coverage","{sample}.tss_histogram.csv"),
    params:
        noise_upper = ( config["filtering"]["tss_slop"] * 2 ) - config["filtering"]["noise_lower"],
        double_slop = ( config["filtering"]["tss_slop"] * 2 ),
        genome_size = config["refs"]["genome_size_bp"],
        tss_slop = config["filtering"]["tss_slop"],
        noise_lower = config["filtering"]["noise_lower"],
    resources:
        mem_mb=16000,
        runtime = 30,
    threads: config["resources"].get("threads", 2)
    conda:
        "../envs/pybedtools.yaml",
    log:
        "logs/rules/tss_coverage_{sample}.log"
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
        bam = os.path.join(result_path, "important_processed", "bam", "{sample}.filtered.bam"),
        bai = os.path.join(result_path, "important_processed", "bam", "{sample}.filtered.bam.bai"),
        peak_file = os.path.join(result_path, "important_processed", "peaks", "{sample}_peaks.narrowPeak"),
        tss_file = os.path.join(result_path, "refs", "tss.bed"),
        excl_regs_file = config["refs"]["blacklist"],
        autosom_ref_file = os.path.join(result_path, "refs", "autosomes.txt")
    output:
        json = os.path.join(result_path, "report", "ataqv", "{sample}.ataqv.json")
    params:
        organism = config["project"].get("genome", "hg38"),
        mito_name = config["refs"].get("mito_name", "chrM"),
        args = config.get("ataqv_args", "--ignore-read-groups"),
        prefix = "{sample}"
    log:
        "logs/rules/ataqv/{sample}.log"
    resources:
        runtime = 100,
    conda:
        "../envs/ataqv.yaml"
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
        jsons = expand(os.path.join(result_path, "report", "ataqv", "{sample}.ataqv.json"), sample=samples.keys())
    output:
        html = directory(os.path.join(result_path, "report", "ataqv_report"))
    params:
        args = config.get("mkarv_args", "")
    log:
        "logs/rules/mkarv/mkarv.log"
    conda:
        "../envs/ataqv.yaml"
    resources:
        mem_mb=64000,
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


rule bamReproducibility:
    """
    Global reproducibility / correlation heatmap using deeptools:
    - multiBamSummary bins
    - plotCorrelation heatmap (Spearman)

    This is a Snakemake translation of the provided SLURM script.
    """
    input:
        bams = lambda w: [
            os.path.join(result_path, "important_processed", "bam", f"{s}.filtered.bam")
            for s in get_reproducibility_sample(w.sample_rep)
        ]
    output:
        npz = temp(os.path.join(result_path, "report", "bamReproducibility", "{sample_rep}_global_rep.npz")),
        heatmap = os.path.join(result_path, "report", "bamReproducibility", "{sample_rep}_global_rep_heatmap.pdf"),
        matrix = os.path.join(result_path, "report", "bamReproducibility", "{sample_rep}_global_rep_cor.txt"),
    params:
        bin_size = 1000,
        prefix = lambda w: w.sample_rep,
    log:
        "logs/rules/bamReproducibility/{sample_rep}_global_rep.log"
    conda:
        "../envs/dtools.yaml"
    resources:
        mem_mb = 28000,
        runtime = 480,
    threads: 20
    shell:
        """
        mkdir -p $(dirname {output.npz})
        mkdir -p $(dirname {log})

        echo "Starting multiBamSummary..." >&2
        multiBamSummary bins \
          --bamfiles {input.bams} \
          --binSize {params.bin_size} \
          --numberOfProcessors {threads} \
          --outFileName {output.npz} 2> {log}

        echo "Starting plotCorrelation..." >&2
        plotCorrelation \
            --corData {output.npz} \
            --corMethod spearman \
            --whatToPlot heatmap \
            --plotTitle "Spearman Correlation of {params.prefix}" \
            --plotFile {output.heatmap} \
            --outFileCorMatrix {output.matrix} 2>> {log}
        """
