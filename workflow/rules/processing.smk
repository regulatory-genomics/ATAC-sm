# alignment with botwtie2 & samtools
rule align:
    input:
        fasta_fwd = os.path.join(result_path, "trimmed", "{sample}_1.fq.gz"),
        fasta_rev = os.path.join(result_path, "trimmed", "{sample}_2.fq.gz"),
        bowtie2_index = os.path.dirname(config["bowtie2_index"]),
        adapter_fasta = config["adapter_fasta"] if config["adapter_fasta"]!="" else [],
        whitelisted_regions = config["whitelisted_regions"],
    output:
        bam = temp(os.path.join(result_path,"results","{sample}","mapped", "{sample}.bam")),
        output_bai =  temp(os.path.join(result_path,"results","{sample}","mapped", "{sample}.bam.bai")),
        filtered_bam = os.path.join(result_path,"results","{sample}","mapped", "{sample}.filtered.bam"),
        filtered_bai = os.path.join(result_path,"results","{sample}","mapped", "{sample}.filtered.bam.bai"),
        bowtie_log = os.path.join(result_path, 'results', "{sample}", 'mapped', '{sample}.txt'),
        bowtie_met = os.path.join(result_path, 'results', "{sample}", 'mapped', '{sample}.bowtie2.met'),
        samblaster_log = os.path.join(result_path, 'results', "{sample}", 'mapped', '{sample}.samblaster.log'),
        samtools_log = os.path.join(result_path, 'results', "{sample}", 'mapped', '{sample}.samtools.log'),
        samtools_flagstat_log = os.path.join(result_path, 'results', "{sample}", 'mapped', '{sample}.samtools_flagstat.log'),
        stats = os.path.join(result_path, 'results', "{sample}", '{sample}.align.stats.tsv'),
    params:
        bowtie2_input = lambda w, input: f"-1 {input.fasta_fwd} -2 {input.fasta_rev}" if samples[f"{w.sample}"]["read_type"] == "paired" else f"-U {input.fasta_fwd}",
        filtering = lambda w: "-q 30 -F {flag} -f 2 -L {whitelist}".format(flag=config['SAM_flag'], whitelist=config["whitelisted_regions"]) if samples["{}".format(w.sample)]["read_type"] == "paired" else "-q 30 -F {flag} -L {whitelist}".format(flag=config['SAM_flag'], whitelist=config["whitelisted_regions"]),
        add_mate_tags = lambda w: "--addMateTags" if samples["{}".format(w.sample)]["read_type"] == "paired" else " ",
        adapter_sequence = "-a " + config["adapter_sequence"] if config["adapter_sequence"] !="" else " ",
        adapter_fasta = "--adapter_fasta " + config["adapter_fasta"] if config["adapter_fasta"] !="" else " ",
        sequencing_platform = config["sequencing_platform"],
        mitochondria_name = config["mitochondria_name"],
        bowtie2_index = config["bowtie2_index"], # The basename of the index for the reference genome excluding the file endings e.g., *.1.bt2
    resources:
        mem_mb=config.get("mem", "16000"),
    threads: 4*config.get("threads", 2)
    conda:
        "../envs/bowtie2.yaml",
    log:
        "logs/rules/align_{sample}.log"
    shell:
        """
        result_path=$(dirname {output.stats})
        find $result_path -type f -name '*.bam.tmp.*' -exec rm {{}} +;
        
        RG="--rg-id {wildcards.sample} --rg SM:{wildcards.sample} --rg PL:{params.sequencing_platform}"

        bowtie2 $RG --very-sensitive --no-discordant -p {threads} --maxins 2000 -x {params.bowtie2_index} \
            --met-file "{output.bowtie_met}" {params.bowtie2_input} 2> "{output.bowtie_log}" | \
            samblaster {params.add_mate_tags} 2> "{output.samblaster_log}" | \
            samtools sort -o "{output.bam}" - 2>> "{output.samtools_log}";
        
        samtools index "{output.bam}" 2>> "{output.samtools_log}";
        samtools idxstats "{output.bam}" | awk '{{ sum += $3 + $4; if($1 == "{params.mitochondria_name}") {{ mito_count = $3; }}}}END{{ print "mitochondrial_fraction\t"mito_count/sum }}' > "{output.stats}";
        samtools flagstat "{output.bam}" > "{output.samtools_flagstat_log}";

        samtools view {params.filtering} -o "{output.filtered_bam}" "{output.bam}";
        samtools index "{output.filtered_bam}";
        """
        
# peak calling with MACS2 & samtools and annotation with HOMER
rule peak_calling:
    input:
        bam = os.path.join(result_path,"results","{sample}","mapped", "{sample}.filtered.bam"),
        bai = os.path.join(result_path,"results","{sample}","mapped", "{sample}.filtered.bam.bai"),
        homer_script = os.path.join(HOMER_path,"configureHomer.pl"),
        regulatory_regions = config["regulatory_regions"],
    output:
        peak_calls = os.path.join(result_path,"results","{sample}","peaks","{sample}_peaks.narrowPeak"),
        peak_annot = os.path.join(result_path,"results","{sample}","peaks","{sample}_peaks.narrowPeak.annotated.tsv"),
        peak_annot_log = os.path.join(result_path,"results","{sample}","peaks","{sample}_peaks.narrowPeak.annotated.tsv.log"),
        macs2_xls = os.path.join(result_path,"results","{sample}","peaks","{sample}_peaks.xls"),
        summits_bed = os.path.join(result_path,"results","{sample}","peaks","{sample}_summits.bed"),
        homer_knownResults = os.path.join(result_path,"results","{sample}","homer","knownResults.txt"),
        homer_log = os.path.join(result_path,"results","{sample}","homer","{sample}.homer.log"),
        macs2_log = os.path.join(result_path, 'results', "{sample}", 'peaks', '{sample}.macs2.log'),
        stats = os.path.join(result_path, 'results', "{sample}", '{sample}.peak.stats.tsv'),
    params:
        peaks_dir = os.path.join(result_path,"results","{sample}","peaks"),
        homer_dir = os.path.join(result_path,"results","{sample}","homer"),
        homer_bin = os.path.join(HOMER_path,"bin"),
        formating = lambda w: '--format BAMPE' if samples["{}".format(w.sample)]["read_type"] == "paired" else '--format BAM',
        genome_size = config["genome_size"],
        genome = config["genome"],
        keep_dup = config['macs2_keep_dup'],
    resources:
        mem_mb=config.get("mem", "16000"),
    threads: 4*config.get("threads", 2)
    conda:
        "../envs/macs2_homer.yaml",
    log:
        "logs/rules/peak_calling_{sample}.log"
    shell:
        """
        export PATH="{params.homer_bin}:$PATH";
        
        macs2 callpeak -t {input.bam} {params.formating} \
            --nomodel --keep-dup {params.keep_dup} --extsize 147 -g {params.genome_size} \
            -n {wildcards.sample} \
            --outdir {params.peaks_dir} > "{output.macs2_log}" 2>&1;
        
        {params.homer_bin}/annotatePeaks.pl {output.peak_calls} {params.genome} > {output.peak_annot} 2> {output.peak_annot_log};
        
        {params.homer_bin}/findMotifsGenome.pl "{output.summits_bed}" {params.genome} {params.homer_dir} -size 200 -mask > "{output.homer_log}" 2>&1

        cat {output.peak_calls} | wc -l | awk '{{print "peaks\t" $1}}' >> "{output.stats}"
        
        TOTAL_READS=`samtools idxstats {input.bam} | awk '{{sum += $3}}END{{print sum}}'`;
        
        samtools view -c -L {output.peak_calls} {input.bam} | awk -v total=$TOTAL_READS '{{print "frip\t" $1/total}}' >> "{output.stats}";

        samtools view -c -L {input.regulatory_regions} {input.bam} | awk -v total=$TOTAL_READS '{{print "regulatory_fraction\t" $1/total}}' >> "{output.stats}";
        
        if [ ! -f {output.homer_knownResults} ]; then
            touch {output.homer_knownResults}
        fi
        """
        
rule aggregate_stats:
    input:
        align_stats = os.path.join(result_path, 'results', "{sample}", '{sample}.align.stats.tsv'),
        peak_stats = os.path.join(result_path, 'results', "{sample}", '{sample}.peak.stats.tsv'),
    output:
        os.path.join(result_path, 'results', "{sample}", '{sample}.stats.tsv'),
    resources:
        mem_mb=config.get("mem", "1000"),
    threads: config.get("threads", 2)
    log:
        "logs/rules/aggregate_stats_{sample}.log"
    shell:
        """
        cat {input.align_stats} {input.peak_stats} > {output}
        """
