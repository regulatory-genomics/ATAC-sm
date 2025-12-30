# ============================================================================
# Configuration & Constants
# ============================================================================
from pathlib import Path

prealign_enabled = config["alignment"].get("prealign", {}).get("enabled", True)
prealignments = config["alignment"].get("prealign", {}).get("indices", []) or []
has_prealignments = prealign_enabled and len(prealignments) > 0
ALIGNER_TOOL = config["alignment"].get("tool", "bowtie2")


# ============================================================================
# Prealignment Rule
# ============================================================================

if has_prealignments:
    # Validate prealignment configuration
    for entry in prealignments:
        if isinstance(entry, dict):
            if "path" not in entry:
                raise ValueError("Each prealignment must define a 'path' key.")
        elif not isinstance(entry, str):
            raise ValueError("Prealignment entries must be dictionaries or strings of the form 'name=index'.")

    rule prealign_reads:
        input:
            fq1 = lambda w: get_trimmed_fastq_paths(w.sample_run)[0],
            fq2 = lambda w: get_trimmed_fastq_paths(w.sample_run)[1] if annot.loc[w.sample_run, "read_type"] == "paired" else [],
        wildcard_constraints:
            sample_run="|".join(annot.index.tolist())
        output:
            fq1 = temp(str(get_output_dir("middle_files/prealigned") / "{sample_run}_prealigned_1.fq.gz")),
            fq2 = temp(str(get_output_dir("middle_files/prealigned") / "{sample_run}_prealigned_2.fq.gz")),
            stats = str(get_output_dir("report/prealigned") / "{sample_run}.prealign.stats.tsv"),
        params:
            prealignments = prealignments,
            aligner = ALIGNER_TOOL,
            prealign_dir = str(get_output_dir("middle_files/prealigned")),
            is_paired = lambda w: annot.loc[w.sample_run, "read_type"] == "paired",
        resources:
            mem_mb = _bwa_mem_mb,
            runtime = 300,
        threads: 5 * config["resources"].get("threads", 2)
        conda:
            "../envs/bwa.yaml" if ALIGNER_TOOL == "bwa-mem2" else "../envs/bowtie2.yaml"
        log:
            "logs/rules/prealign_{sample_run}.log"
        script:
            "../scripts/prealign.py"

# ============================================================================
# Alignment Rules
# ============================================================================

if ALIGNER_TOOL == "bowtie2":
    rule align_bowtie2:
        input:
            fasta_fwd = lambda w: get_reads(w, 0),
            fasta_rev = lambda w: get_reads(w, 1),
            bowtie2_index = os.path.dirname(config["alignment"]["bowtie2"]["index"]),
            adapter_fasta = config["adapters"]["fasta"] if config["adapters"]["fasta"] != "" else [],
            whitelisted_regions = config["refs"]["whitelist"],
        wildcard_constraints:
            sample="|".join(samples.keys())
        output:
            bam = str(get_output_dir("important_processed/bam") / "{sample}.filtered.bam"),
            bai = str(get_output_dir("important_processed/bam") / "{sample}.filtered.bam.bai"),
            bowtie_log = str(get_output_dir("logs/align") / "{sample}.bowtie2.log"),
            bowtie_met = str(get_output_dir("logs/align") / "{sample}.bowtie2.met"),
            samblaster_log = str(get_output_dir("logs/align") / "{sample}.samblaster.log"),
        params:
            sample_name = "{sample}",
            # Wrapped in lambda for safety
            bowtie2_input = lambda w, input: get_bowtie2_input_string(w, input),
            add_mate_tags = get_add_mate_tags,
            adapter_sequence = "-a " + config["adapters"]["sequence"] if config["adapters"]["sequence"] != "" else "",
            adapter_fasta = "--adapter_fasta " + config["adapters"]["fasta"] if config["adapters"]["fasta"] != "" else "",
            sequencing_platform = config["alignment"]["sequencing_platform"],
            bowtie2_index = config["alignment"]["bowtie2"]["index"],
            bowtie2_local_mode = "--local" if config["alignment"].get("local_mode", False) else "",
            filtering_flags = get_filtering_flags,
        resources:
            mem_mb = config["resources"].get("mem_mb", 16000),
            runtime = 1000,
        threads: 5 * config["resources"].get("threads", 2)
        conda:
            "../envs/bowtie2.yaml"
        log:
            "logs/rules/align_{sample}.log"
        shell:
            """
            set -euo pipefail
            
            mkdir -p $(dirname {output.bam}) $(dirname {output.bowtie_log})
            result_path=$(dirname {output.bam})
            find $result_path -type f -name '{wildcards.sample}.filtered.bam.tmp.*' -delete 2>/dev/null || true
            rm -f "{output.bowtie_log}" "{output.bowtie_met}" "{output.samblaster_log}" 2>/dev/null || true
            
            RG="--rg-id {wildcards.sample} --rg SM:{params.sample_name} --rg PL:{params.sequencing_platform}"
            
            bowtie2 $RG --very-sensitive --no-discordant -p {threads} --maxins 2000 \
                -x {params.bowtie2_index} {params.bowtie2_local_mode} \
                --met-file "{output.bowtie_met}" {params.bowtie2_input} 2> "{output.bowtie_log}" | \
            samblaster {params.add_mate_tags} 2> "{output.samblaster_log}" | \
            samtools view {params.filtering_flags} -bhS - 2>> "{output.bowtie_log}" | \
            samtools sort -o "{output.bam}" - 2>> "{output.bowtie_log}"
            
            samtools index "{output.bam}" 2>> "{output.bowtie_log}"
            """

elif ALIGNER_TOOL == "bwa-mem2":
    rule align_bwa_mem:
        input:
            fasta_fwd = lambda w: get_reads(w, 0),
            fasta_rev = lambda w: get_reads(w, 1),
            index = get_bwa_index_input,
            whitelisted_regions = config["refs"]["whitelist"],
        wildcard_constraints:
            sample="|".join(samples.keys())
        output:
            bam = str(get_output_dir("important_processed/bam") / "{sample}.filtered.bam"),
            bai = str(get_output_dir("important_processed/bam") / "{sample}.filtered.bam.bai"),
            bwa_log = str(get_output_dir("logs/align/{sample}") / "{sample}.bwa.log"),
            samblaster_log = str(get_output_dir("logs/align/{sample}") / "{sample}.samblaster.log"),
        params:
            sample_name = "{sample}",
            # Wrapped in lambda for safety
            bwa_input = lambda w, input: get_bwa_input_string(w, input),
            add_mate_tags = get_add_mate_tags,
            sequencing_platform = config["alignment"]["sequencing_platform"],
            bwa_args = config["alignment"]["bwa"].get("extra_args", ""),
            bwa_min_score_flag = sanitize_bwa_min_score_flag(),
            bwa_m_flag = "-M",
            bwa_index_path = get_bwa_index_path(),
            filtering_flags = get_filtering_flags,
        resources:
            mem_mb = _bwa_mem_mb,
            runtime = 800,
        threads: 5 * config["resources"].get("threads", 2)
        conda:
            "../envs/bwa.yaml"
        log:
            "logs/rules/align_bwa_mem_{sample}.log"
        shell:
            """
            set -euo pipefail
            
            mkdir -p $(dirname {output.bam}) $(dirname {output.bwa_log})
            result_path=$(dirname {output.bam})
            find $result_path -type f -name '{wildcards.sample}.filtered.bam.tmp.*' -delete 2>/dev/null || true
            rm -f "{output.bwa_log}" "{output.samblaster_log}" 2>/dev/null || true
            
            RG="@RG\\tID:{wildcards.sample}\\tSM:{params.sample_name}\\tPL:{params.sequencing_platform}"
            
            bwa-mem2 mem \
                {params.bwa_args} \
                {params.bwa_m_flag} \
                {params.bwa_min_score_flag} \
                -R "$RG" \
                -t {threads} \
                {params.bwa_index_path} \
                {params.bwa_input} 2> {output.bwa_log} | \
            samblaster {params.add_mate_tags} 2> {output.samblaster_log} | \
            samtools view {params.filtering_flags} -bhS - 2>> {output.bwa_log} | \
            samtools sort -o {output.bam} - 2>> {output.bwa_log}
            
            samtools index "{output.bam}" 2>> "{output.bwa_log}"
            """

else:
    raise ValueError(f"Unknown aligner: {ALIGNER_TOOL}. Must be 'bowtie2' or 'bwa-mem2'")

# ============================================================================
# BWA Index Rule
# ============================================================================

rule bwa_mem2_index:
    input:
        fasta = config["refs"]["fasta"]
    output:
        index = multiext(config["refs"]["fasta"], ".amb", ".ann", ".pac", ".bwt.2bit.64", ".0123"),
    log:
        "logs/bwa_mem2_index/bwa_mem2_index.log"
    conda:
        "../envs/bwa.yaml"
    shell:
        """
        bwa-mem2 index {input.fasta} 2> {log}
        """

# ============================================================================
# Alignment Statistics
# ============================================================================

rule samtools_process:
    input:
        bam = str(get_output_dir("important_processed/bam") / "{sample}.filtered.bam"),
        bai = str(get_output_dir("important_processed/bam") / "{sample}.filtered.bam.bai"),
    output:
        samtools_log = str(get_output_dir("logs/align") / "{sample}.samtools.log"),
        samtools_flagstat_log = str(get_output_dir("logs/align") / "{sample}.samtools_flagstat.log"),
        stats = str(get_output_dir("report/align_stats") / "{sample}.align.stats.tsv"),
    params:
        mitochondria_name = config["refs"].get("mito_name", "chrM"),
    resources:
        mem_mb = config["resources"].get("mem_mb", 16000),
        runtime = 90,
    threads: config["resources"].get("threads", 1)
    conda:
        "../envs/bowtie2.yaml"
    wildcard_constraints:
        sample="|".join(samples.keys())
    log:
        "logs/rules/samtools_process_{sample}.log"
    shell:
        """
        set -euo pipefail
        
        mkdir -p $(dirname {output.stats})
        samtools idxstats "{input.bam}" | awk '{{ 
            sum += $3 + $4; 
            if($1 == "{params.mitochondria_name}") {{ mito_count = $3; }}
        }}END{{ 
            print "mitochondrial_fraction\\t"mito_count/sum 
        }}' > "{output.stats}" 2>> "{output.samtools_log}"
        
        samtools flagstat "{input.bam}" > "{output.samtools_flagstat_log}" 2>> "{output.samtools_log}"
        """

# ============================================================================
# BAM to TagAlign Conversion
# ============================================================================

rule bam_to_bed:
    input:
        bam = str(get_output_dir("important_processed/bam") / "{sample}.filtered.bam"),
        bai = str(get_output_dir("important_processed/bam") / "{sample}.filtered.bam.bai"),
    output:
        bed = str(get_output_dir("middle_files/bed") / "{sample}.tagAlign.gz"),
    params:
        script = os.path.join(workflow.basedir, "scripts", "bam_to_tagalign.sh"),
        bed_dir = str(get_output_dir("middle_files/bed")),
        is_paired = lambda w: samples[w.sample].get("read_type", "single") == "paired",
        disable_tn5_shift = config.get("peaks", {}).get("disable_tn5_shift", False),
    resources:
        mem_mb = config["resources"].get("mem_mb", 16000),
        runtime = 300,
    threads: config["resources"].get("threads", 2)
    conda:
        "../envs/pybedtools.yaml"
    log:
        "logs/rules/bam_to_bed_{sample}.log"
    wildcard_constraints:
        sample="|".join(samples.keys())
    shell:
        """
        {params.script} \\
            {input.bam} \\
            {output.bed} \\
            {wildcards.sample} \\
            {params.is_paired} \\
            {params.disable_tn5_shift} \\
            {threads} \\
            {params.bed_dir} \\
            {log}
        """

# ============================================================================
# Peak Calling
# ============================================================================

rule peak_calling:
    input:
        bed = str(get_output_dir("middle_files/bed") / "{sample}.tagAlign.gz"),
    output:
        peak_calls = str(get_output_dir("important_processed/peaks") / "{sample}_peaks.narrowPeak"),
        macs2_xls = str(get_output_dir("important_processed/peaks") / "{sample}_peaks.xls"),
        summits_bed = str(get_output_dir("important_processed/peaks") / "{sample}_summits.bed"),
        macs2_log = str(get_output_dir("important_processed/peaks") / "{sample}.macs2.log"),
    params:
        peaks_dir = str(get_output_dir("important_processed/peaks")),
        genome_size = config["refs"]["genome_size_bp"],
        keep_dup = config["peaks"]["macs2_keep_dup"],
        macs2_shift = config["peaks"].get("macs2_shift", -75),
        macs2_extsize = config["peaks"].get("macs2_extsize", 150),
        macs2_format = config["peaks"].get("macs2_format", "BED"),
        pval = config["peaks"].get("macs2_pval", 1e-3),
    resources:
        mem_mb = config["resources"].get("mem_mb", 16000),
        runtime = 60,
    threads: config["resources"].get("threads", 2)
    conda:
        "../envs/macs2_homer.yaml"
    log:
        "logs/rules/peak_calling_{sample}.log"
    wildcard_constraints:
        sample="|".join(samples.keys())
    shell:
        """
        set -euo pipefail
        
        mkdir -p {params.peaks_dir}
        
        macs2 callpeak -t {input.bed} -f {params.macs2_format} \
            --nomodel --keep-dup {params.keep_dup} \
            --shift {params.macs2_shift} --extsize {params.macs2_extsize} \
            -g {params.genome_size} \
            -n {wildcards.sample} \
            -p {params.pval} \
            --outdir {params.peaks_dir} > "{output.macs2_log}" 2>&1
        
        if [ ! -f {output.peak_calls} ] || [ ! -s {output.peak_calls} ]; then
            touch {output.peak_calls}
            touch {output.summits_bed}
            echo "Warning: No peaks called for {wildcards.sample}" >> "{output.macs2_log}"
        fi
        """

# ============================================================================
# Peak Annotation
# ============================================================================

rule peak_annotation:
    input:
        peak_calls = str(get_output_dir("important_processed/peaks") / "{sample}_peaks.narrowPeak"),
        summits_bed = str(get_output_dir("important_processed/peaks") / "{sample}_summits.bed"),
        bam = str(get_output_dir("important_processed/bam") / "{sample}.filtered.bam"),
        bai = str(get_output_dir("important_processed/bam") / "{sample}.filtered.bam.bai"),
        homer_script = os.path.join(HOMER_path, "configureHomer.pl"),
        regulatory_regions = config["refs"]["regulatory_regions"],
    output:
        peak_annot = str(get_output_dir("important_processed/peaks") / "{sample}_peaks.narrowPeak.annotated.tsv"),
        peak_annot_log = str(get_output_dir("important_processed/peaks") / "{sample}_peaks.narrowPeak.annotated.tsv.log"),
        homer_knownResults = str(get_output_dir("downstream_res/homer/{sample}") / "knownResults.txt"),
        homer_log = str(get_output_dir("downstream_res/homer") / "{sample}.homer.log"),
        stats = str(get_output_dir("report/peaks_stats") / "{sample}.stats.tsv"),
    params:
        peaks_dir = str(get_output_dir("important_processed/peaks")),
        homer_dir = str(get_output_dir("downstream_res/homer/{sample}")),
        homer_bin = os.path.join(HOMER_path, "bin"),
        genome = config["project"]["genome"],
    resources:
        mem_mb = config["resources"].get("mem_mb", 16000),
        runtime = 240,
    threads: config["resources"].get("threads", 2)
    conda:
        "../envs/macs2_homer.yaml"
    log:
        "logs/rules/peak_annotation_{sample}.log"
    wildcard_constraints:
        sample="|".join(samples.keys())
    shell:
        """
        set -euo pipefail
        
        export PATH="{params.homer_bin}:$PATH"
        mkdir -p {params.homer_dir}
        
        # Initialize output files
        touch {output.peak_annot} {output.peak_annot_log} {output.homer_log} {output.homer_knownResults} {output.stats}
        
        if [ -s "{input.peak_calls}" ]; then
            # Annotate peaks
            {params.homer_bin}/annotatePeaks.pl {input.peak_calls} {params.genome} > {output.peak_annot} 2> {output.peak_annot_log} || \
                echo "HOMER annotation completed with warnings" >> {output.peak_annot_log}
            
            # Motif finding
            if [ -s "{input.summits_bed}" ]; then
                {params.homer_bin}/findMotifsGenome.pl "{input.summits_bed}" {params.genome} {params.homer_dir} \
                    -size 200 -mask > "{output.homer_log}" 2>&1 || \
                    echo "HOMER motif finding completed with warnings" >> {output.homer_log}
            else
                echo "No summits file available for motif finding" > {output.homer_log}
            fi
            
            # Calculate statistics
            PEAK_COUNT=$(wc -l < {input.peak_calls} || echo 0)
            echo -e "peaks\\t$PEAK_COUNT" > {output.stats}
            
            TOTAL_READS=$(samtools idxstats {input.bam} 2>/dev/null | awk '{{sum += $3}}END{{print sum+0}}' || echo "0")
            
            if [ "$TOTAL_READS" -gt 0 ] 2>/dev/null; then
                FRIP=$(samtools view -c -L {input.peak_calls} {input.bam} 2>/dev/null | \
                    awk -v total=$TOTAL_READS '{{if(NF>0) print $1/total; else print 0}}' || echo "0")
                echo -e "frip\\t$FRIP" >> {output.stats}
                
                REGULATORY_FRAC=$(samtools view -c -L {input.regulatory_regions} {input.bam} 2>/dev/null | \
                    awk -v total=$TOTAL_READS '{{if(NF>0) print $1/total; else print 0}}' || echo "0")
                echo -e "regulatory_fraction\\t$REGULATORY_FRAC" >> {output.stats}
            else
                echo -e "frip\\t0" >> {output.stats}
                echo -e "regulatory_fraction\\t0" >> {output.stats}
            fi
        else
            echo "No peaks found for {wildcards.sample}, creating empty annotation files" > {output.peak_annot_log}
            echo "No peaks found for {wildcards.sample}, skipping HOMER motif finding" > {output.homer_log}
            echo -e "peaks\\t0" > {output.stats}
            echo -e "frip\\t0" >> {output.stats}
            echo -e "regulatory_fraction\\t0" >> {output.stats}
        fi
        
        if [ ! -f {output.homer_knownResults} ]; then
            touch {output.homer_knownResults}
        fi
        """

# ============================================================================
# BigWig Track Generation
# ============================================================================

rule tracks:
    input:
        bam = str(get_output_dir("important_processed/bam") / "{sample}.filtered.bam"),
        bai = str(get_output_dir("important_processed/bam") / "{sample}.filtered.bam.bai"),
    output:
        bw = str(get_output_dir("important_processed/tracks") / "{sample}.bw"),
    conda:
        "../envs/dtools.yaml"
    resources:
        mem_mb = config["resources"].get("mem_mb", 16000),
        runtime = 120,
    threads: 8
    wildcard_constraints:
        sample="|".join(samples.keys())
    log:
        "logs/rules/tracks_{sample}.log"
    shell:
        """
        set -euo pipefail
        
        mkdir -p $(dirname {output.bw})
        bamCoverage -b {input.bam} -o {output.bw} \
            --binSize 10 --smoothLength 50 --normalizeUsing CPM \
            -p {threads} 2> {log}
        """
