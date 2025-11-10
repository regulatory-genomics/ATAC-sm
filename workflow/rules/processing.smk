prealign_enabled = config.get("prealign_enabled", True)
prealignments = config.get("prealignments", []) or []
has_prealignments = prealign_enabled and len(prealignments) > 0


def _trimmed_fastq_paths(sample_run):
    if annot.loc[sample_run, "read_type"] == "paired":
        return (
            os.path.join(result_path, "trimmed", f"{sample_run}_1.fq.gz"),
            os.path.join(result_path, "trimmed", f"{sample_run}_2.fq.gz"),
        )
    return (
        os.path.join(result_path, "trimmed", f"{sample_run}.fq.gz"),
        None,
    )


def _prealigned_fastq_paths(sample_run):
    base_dir = os.path.join(result_path, "results", sample_run, "prealign")
    return (
        os.path.join(base_dir, f"{sample_run}_prealigned_1.fq.gz"),
        os.path.join(base_dir, f"{sample_run}_prealigned_2.fq.gz"),
    )


if has_prealignments:
    for entry in prealignments:
        if isinstance(entry, dict):
            if "index" not in entry:
                raise ValueError("Each prealignment must define an 'index' key.")
        elif not isinstance(entry, str):
            raise ValueError("Prealignment entries must be dictionaries or strings of the form 'name=index'.")

    rule prealign_reads:
        input:
            fq1 = lambda w: _trimmed_fastq_paths(w.sample_run)[0],
            fq2 = lambda w: _trimmed_fastq_paths(w.sample_run)[1] if annot.loc[w.sample_run, "read_type"] == "paired" else [],
        wildcard_constraints:
            sample_run="|".join(annot.index.tolist())
        output:
            fq1 = temp(os.path.join(result_path, "results", "{sample_run}", "prealign", "{sample_run}_prealigned_1.fq.gz")),
            fq2 = temp(os.path.join(result_path, "results", "{sample_run}", "prealign", "{sample_run}_prealigned_2.fq.gz")),
            stats = os.path.join(result_path, "results", "{sample_run}", "prealign", "{sample_run}.prealign.stats.tsv"),
        params:
            prealignments = prealignments,
            aligner = config.get("aligner", "bowtie2"),
            prealign_dir = lambda w: os.path.join(result_path, "results", w.sample_run, "prealign"),
            is_paired = lambda w: annot.loc[w.sample_run, "read_type"] == "paired",
        resources:
            mem_mb = config.get("mem", "16000"),
        threads: config.get("threads", 2)
        conda:
            "../envs/bwa.yaml" if config.get("aligner", "bowtie2") == "bwa-mem2" else "../envs/bowtie2.yaml"
        log:
            "logs/rules/prealign_{sample_run}.log"
        run:
            import os
            import shutil
            import gzip
            from pathlib import Path
            from shlex import quote
            from snakemake.shell import shell

            prealign_list = params.prealignments
            log_file = str(log)
            os.makedirs(os.path.dirname(log_file), exist_ok=True)
            with open(log_file, "w"):
                pass

            prealign_dir = params.prealign_dir
            os.makedirs(prealign_dir, exist_ok=True)

            paired = params.is_paired
            original_fq1 = str(input.fq1)
            original_fq2 = str(input.fq2) if paired and input.fq2 else None
            current_fq1 = original_fq1
            current_fq2 = original_fq2
            aligner = params.aligner.lower()

            stats_path = str(output.stats)
            stats_lines = []

            def count_reads(path):
                if not path or not os.path.exists(path):
                    return 0
                opener = gzip.open if path.endswith(".gz") else open
                total_lines = 0
                with opener(path, "rt") as handle:
                    for _ in handle:
                        total_lines += 1
                return total_lines // 4

            output_fq1 = str(output.fq1)
            output_fq2 = str(output.fq2)

            if not prealign_list:
                shutil.copyfile(current_fq1, output_fq1)
                if paired and current_fq2:
                    shutil.copyfile(current_fq2, output_fq2)
                else:
                    Path(output_fq2).parent.mkdir(parents=True, exist_ok=True)
                    Path(output_fq2).touch(exist_ok=True)
                total_reads = count_reads(current_fq1)
                with open(stats_path, "w") as stats_fh:
                    stats_fh.write("prealignment\treads_before\treads_after\treads_filtered\tpercent_filtered\n")
                    stats_fh.write(f"total\t{total_reads}\t{total_reads}\t0\t0.000000\n")
                return

            def _parse_entry(entry):
                if isinstance(entry, dict):
                    idx = entry.get("index")
                    if idx is None:
                        raise ValueError("Prealignment dictionary entries require an 'index'.")
                    name = entry.get("name") or Path(idx).stem
                    return name, idx
                entry = str(entry)
                if "=" in entry:
                    name, idx = entry.split("=", 1)
                    return name.strip(), idx.strip()
                return Path(entry).stem, entry

            def _cleanup(path):
                if not path or not isinstance(path, str):
                    return
                keep = {original_fq1, output_fq1}
                if paired and original_fq2:
                    keep.add(original_fq2)
                if paired and output_fq2:
                    keep.add(output_fq2)
                if path.startswith(prealign_dir) and os.path.exists(path) and path not in keep:
                    try:
                        os.remove(path)
                    except:
                        pass

            for entry in prealign_list:
                name, index = _parse_entry(entry)
                stage_prefix = os.path.join(prealign_dir, f"{wildcards.sample_run}.{name}.unmapped")
                prev_fq1 = current_fq1
                prev_fq2 = current_fq2
                reads_before = count_reads(prev_fq1)

                if aligner == "bowtie2":
                    if paired:
                        cmd = (
                            "set -euo pipefail; "
                            f"bowtie2 -x {quote(index)} -1 {quote(current_fq1)} -2 {quote(current_fq2)} "
                            f"-p {threads} --un-conc-gz {quote(stage_prefix)} -S /dev/null "
                            f"2>> {quote(log_file)}"
                        )
                        shell(cmd)
                        current_fq1 = f"{stage_prefix}.1.gz"
                        current_fq2 = f"{stage_prefix}.2.gz"
                    else:
                        single_target = f"{stage_prefix}.gz"
                        cmd = (
                            "set -euo pipefail; "
                            f"bowtie2 -x {quote(index)} -U {quote(current_fq1)} "
                            f"-p {threads} --un-gz {quote(single_target)} -S /dev/null "
                            f"2>> {quote(log_file)}"
                        )
                        shell(cmd)
                        current_fq1 = single_target
                        current_fq2 = None
                else:
                    # bwa-mem2
                    if paired:
                        fq1_tmp = f"{stage_prefix}_1.fq"
                        fq2_tmp = f"{stage_prefix}_2.fq"
                        cmd = (
                            "set -euo pipefail; "
                            f"bwa-mem2 mem -t {threads} {quote(index)} {quote(current_fq1)} {quote(current_fq2)} "
                            f"2>> {quote(log_file)} | "
                            f"samtools fastq -f 4 -1 {quote(fq1_tmp)} -2 {quote(fq2_tmp)} "
                            f"-0 /dev/null -s /dev/null -n - 2>> {quote(log_file)}"
                        )
                        shell(cmd)
                        shell(f"gzip -f {quote(fq1_tmp)}")
                        shell(f"gzip -f {quote(fq2_tmp)}")
                        current_fq1 = f"{fq1_tmp}.gz"
                        current_fq2 = f"{fq2_tmp}.gz"
                    else:
                        fq_tmp = f"{stage_prefix}.fq"
                        cmd = (
                            "set -euo pipefail; "
                            f"bwa-mem2 mem -t {threads} {quote(index)} {quote(current_fq1)} "
                            f"2>> {quote(log_file)} | "
                            f"samtools fastq -f 4 -0 {quote(fq_tmp)} -s /dev/null -n - "
                            f"2>> {quote(log_file)}"
                        )
                        shell(cmd)
                        shell(f"gzip -f {quote(fq_tmp)}")
                        current_fq1 = f"{fq_tmp}.gz"
                        current_fq2 = None

                reads_after = count_reads(current_fq1)
                filtered_reads = max(reads_before - reads_after, 0)
                percent_filtered = filtered_reads / reads_before if reads_before else 0.0
                stats_lines.append((name, reads_before, reads_after, filtered_reads, percent_filtered))

                if prev_fq1 != original_fq1:
                    _cleanup(prev_fq1)
                if paired and prev_fq2 and prev_fq2 != original_fq2:
                    _cleanup(prev_fq2)

            os.makedirs(os.path.dirname(output_fq1), exist_ok=True)
            os.makedirs(os.path.dirname(output_fq2), exist_ok=True)

            if current_fq1 != output_fq1:
                if os.path.exists(output_fq1):
                    os.remove(output_fq1)
                shutil.move(current_fq1, output_fq1)

            if paired and current_fq2:
                if current_fq2 != output_fq2:
                    if os.path.exists(output_fq2):
                        os.remove(output_fq2)
                    shutil.move(current_fq2, output_fq2)
            else:
                Path(output_fq2).touch(exist_ok=True)

            if stats_lines:
                total_before = stats_lines[0][1]
                total_after = stats_lines[-1][2]
            else:
                total_before = count_reads(original_fq1)
                total_after = count_reads(output_fq1)
            total_filtered = max(total_before - total_after, 0)
            total_percent = total_filtered / total_before if total_before else 0.0

            with open(stats_path, "w") as stats_fh:
                stats_fh.write("prealignment\treads_before\treads_after\treads_filtered\tpercent_filtered\n")
                for row in stats_lines:
                    stats_fh.write(f"{row[0]}\t{row[1]}\t{row[2]}\t{row[3]}\t{row[4]:.6f}\n")
                stats_fh.write(f"total\t{total_before}\t{total_after}\t{total_filtered}\t{total_percent:.6f}\n")


if config.get("aligner", "bowtie2") == "bowtie2":
    # alignment with bowtie2 & samtools (per run)
    rule align_bowtie2:
        input:
            fasta_fwd = lambda w: (_prealigned_fastq_paths(w.sample_run)[0] if has_prealignments else _trimmed_fastq_paths(w.sample_run)[0]),
            fasta_rev = lambda w: (_prealigned_fastq_paths(w.sample_run)[1] if has_prealignments else _trimmed_fastq_paths(w.sample_run)[1]),
            bowtie2_index = os.path.dirname(config["bowtie2_index"]),
            adapter_fasta = config["adapter_fasta"] if config["adapter_fasta"]!="" else [],
            whitelisted_regions = config["whitelisted_regions"],
        wildcard_constraints:
            sample_run="|".join(annot.index.tolist())  # Only match actual sample_run names from annotation
        output:
            bam = temp(os.path.join(result_path,"results","{sample_run}","mapped", "{sample_run}.bam")),
            bowtie_log = os.path.join(result_path, 'results', "{sample_run}", 'mapped', '{sample_run}.txt'),
            bowtie_met = os.path.join(result_path, 'results', "{sample_run}", 'mapped', '{sample_run}.bowtie2.met'),
            samblaster_log = os.path.join(result_path, 'results', "{sample_run}", 'mapped', '{sample_run}.samblaster.log'),
        params:
            sample_name = lambda w: annot.loc[w.sample_run, 'sample_name'],
            bowtie2_input = lambda w, input: f"-1 {input.fasta_fwd} -2 {input.fasta_rev}" if annot.loc[w.sample_run, "read_type"] == "paired" else f"-U {input.fasta_fwd}",
            add_mate_tags = lambda w: "--addMateTags" if annot.loc[w.sample_run, "read_type"] == "paired" else " ",
            adapter_sequence = "-a " + config["adapter_sequence"] if config["adapter_sequence"] !="" else " ",
            adapter_fasta = "--adapter_fasta " + config["adapter_fasta"] if config["adapter_fasta"] !="" else " ",
            sequencing_platform = config["sequencing_platform"],
            bowtie2_index = config["bowtie2_index"], # The basename of the index for the reference genome excluding the file endings e.g., *.1.bt2
            bowtie2_local_mode = "--local" if config.get("local", False) else "",  # Add this param for local mode
        resources:
            mem_mb=config.get("mem", "16000"),
        threads: 4*config.get("threads", 2)
        conda:
            "../envs/bowtie2.yaml",
        log:
            "logs/rules/align_{sample_run}.log"
        shell:
            """
            mkdir -p $(dirname {output.bam})
            result_path=$(dirname {output.bam})
            find $result_path -type f -name '*.bam.tmp.*' -exec rm {{}} +;
            
            RG="--rg-id {wildcards.sample_run} --rg SM:{params.sample_name} --rg PL:{params.sequencing_platform}"

            bowtie2 $RG --very-sensitive --no-discordant -p {threads} --maxins 2000 -x {params.bowtie2_index} \
                {params.bowtie2_local_mode}  \
                --met-file "{output.bowtie_met}" {params.bowtie2_input} 2> "{output.bowtie_log}" | \
                samblaster {params.add_mate_tags} 2> "{output.samblaster_log}" | \
                samtools sort -o "{output.bam}" - 2>> "{output.bowtie_log}";
            """

elif config.get("aligner", "bowtie2") == "bwa-mem2":
    def _sanitize_bwa_min_score_flag():
        value = config.get("bwa_min_score")
        if value is None:
            return ""
        if isinstance(value, str):
            stripped = value.strip()
            if stripped.lower() in {"", "null", "none"}:
                return ""
            value_str = stripped
        else:
            value_str = value
        return f"-T {value_str}"

    # Using bwa-mem2 only

    rule align_bwa_mem:
        input:
            fasta_fwd = lambda w: (_prealigned_fastq_paths(w.sample_run)[0] if has_prealignments else _trimmed_fastq_paths(w.sample_run)[0]),
            fasta_rev = lambda w: (_prealigned_fastq_paths(w.sample_run)[1] if has_prealignments else _trimmed_fastq_paths(w.sample_run)[1]),
            # If bwa_mem2_path is set, use BWA-MEM2 index directory (pre-indexed)
            # If bwa_mem2_path is not set, require bwa_index rule to create index files
            index = lambda w: multiext(config["genome_fasta"], ".amb", ".ann", ".bwt", ".pac", ".sa", ".0123", ".alt") if not config.get("bwa_mem2_path") or config.get("bwa_mem2_path") == "null" or config.get("bwa_mem2_path") == "" else [],
            whitelisted_regions = config["whitelisted_regions"],
        wildcard_constraints:
            sample_run="|".join(annot.index.tolist())  # Only match actual sample_run names from annotation
        output:
            bam = temp(os.path.join(result_path,"results","{sample_run}","mapped", "{sample_run}.bam")),
            bwa_log = os.path.join(result_path, 'results', "{sample_run}", 'mapped', '{sample_run}.bwa.log'),
            samblaster_log = os.path.join(result_path, 'results', "{sample_run}", 'mapped', '{sample_run}.samblaster.log'),
        params:
            sample_name = lambda w: annot.loc[w.sample_run, 'sample_name'],
            bwa_input = lambda w, input: f"{input.fasta_fwd} {input.fasta_rev}" if annot.loc[w.sample_run, "read_type"] == "paired" else f"{input.fasta_fwd}",
            add_mate_tags = lambda w: "--addMateTags" if annot.loc[w.sample_run, "read_type"] == "paired" else " ",
            sequencing_platform = config["sequencing_platform"],
            bwa_args = config.get("bwa_args", ""),
            bwa_min_score_flag = _sanitize_bwa_min_score_flag(),
            bwa_m_flag = "-M",  # Mark shorter split hits as secondary
            # Use bwa-mem2 index path if set, otherwise use genome_fasta (for bwa)
            bwa_index_path = lambda w: (config.get("bwa_mem2_path") if config.get("bwa_mem2_path") not in ("", "null", None) else config["genome_fasta"]),
        resources:
            mem_mb=config.get("mem", "64000"),
        threads: 4*config.get("threads", 2)
        conda:
            "../envs/bwa.yaml",
        log:
            "logs/rules/align_bwa_mem_{sample_run}.log"
        shell:
            """
            mkdir -p $(dirname {output.bam})
            result_path=$(dirname {output.bam})
            find $result_path -type f -name '*.bam.tmp.*' -exec rm {{}} +;
            
            RG="@RG\\tID:{wildcards.sample_run}\\tSM:{params.sample_name}\\tPL:{params.sequencing_platform}"

            BWA_INDEX_PATH="{params.bwa_index_path}"

            bwa-mem2 mem \
                {params.bwa_args} \
                {params.bwa_m_flag} \
                {params.bwa_min_score_flag} \
                -R "$RG" \
                -t {threads} \
                $BWA_INDEX_PATH \
                {params.bwa_input} 2> {output.bwa_log} | \
                samblaster {params.add_mate_tags} 2> {output.samblaster_log} | \
                samtools view -bhS -F 0x0100 -O BAM - 2>> {output.bwa_log} | \
                samtools sort -o {output.bam} - 2>> {output.bwa_log};
            """

else:
    raise ValueError(f"Unknown aligner: {config.get('aligner', 'bowtie2')}. Must be 'bowtie2' or 'bwa-mem2'")

rule bwa_mem2_index:
    input:
        fasta = config["genome_fasta"]
    output:
        index = multiext(config["genome_fasta"], ".amb", ".ann", ".pac", ".bwt.2bit.64", ".0123"),
    log:
        "logs/bwa_mem2_index/bwa_mem2_index.log"
    conda:
        "../envs/bwa.yaml"
    shell:
        """
        # Create BWA index
        bwa-mem2 index {input.fasta} 2> {log}
        """

rule samtools_process:
    input:
        bam = os.path.join(result_path,"results","{sample_run}","mapped", "{sample_run}.bam"),
    output:
        filtered_bam = temp(os.path.join(result_path,"results","{sample_run}","mapped", "{sample_run}.filtered.bam")),
        filtered_bai = temp(os.path.join(result_path,"results","{sample_run}","mapped", "{sample_run}.filtered.bam.bai")),
        samtools_log = os.path.join(result_path, 'results', "{sample_run}", 'mapped', '{sample_run}.samtools.log'),
        samtools_flagstat_log = os.path.join(result_path, 'results', "{sample_run}", 'mapped', '{sample_run}.samtools_flagstat.log'),
        stats = os.path.join(result_path, 'results', "{sample_run}", '{sample_run}.align.stats.tsv'),
    params:
        filtering = lambda w: "-q 30 -F {flag} -f 2 -L {whitelist}".format(flag=config['SAM_flag'], whitelist=config["whitelisted_regions"]) if annot.loc[w.sample_run, "read_type"] == "paired" else "-q 30 -F {flag} -L {whitelist}".format(flag=config['SAM_flag'], whitelist=config["whitelisted_regions"]),
        mitochondria_name = config["mitochondria_name"],
    resources:
        mem_mb=config.get("mem", "16000"),
    threads: config.get("threads", 1)
    conda:
        "../envs/bowtie2.yaml",
    shell:
        '''
            samtools index "{input.bam}" 2>> "{output.samtools_log}";
            samtools idxstats "{input.bam}" | awk '{{ sum += $3 + $4; if($1 == "{params.mitochondria_name}") {{ mito_count = $3; }}}}END{{ print "mitochondrial_fraction\t"mito_count/sum }}' > "{output.stats}";
            samtools flagstat "{input.bam}" > "{output.samtools_flagstat_log}";

            samtools view {params.filtering} -o "{output.filtered_bam}" "{input.bam}";
            samtools index "{output.filtered_bam}";
        '''

# Merge BAM files for samples with multiple runs
# Note: This rule only matches actual sample names (not sample_run names)
def get_merge_bam_inputs(wildcards):
    """Get inputs for merge_bam, validating that wildcards.sample is a valid sample name"""
    if wildcards.sample not in samples:
        # Return a non-existent file to prevent rule matching for invalid sample names
        # This file will never exist, so the rule won't match
        return [os.path.join(result_path, "__INVALID_SAMPLE__", f"{wildcards.sample}.bam")]
    return [os.path.join(result_path, "results", sr, "mapped", f"{sr}.filtered.bam") 
            for sr in get_runs_for_sample(wildcards.sample)]

rule merge_bam:
    input:
        get_merge_bam_inputs,
    wildcard_constraints:
        sample="|".join(samples.keys())  # Only match actual sample names (not sample_run names)
    output:
        bam = os.path.join(result_path,"bam","{sample}", "{sample}.filtered.bam"),
        bai = os.path.join(result_path,"bam","{sample}", "{sample}.filtered.bam.bai"),
    params:
        sample_name = "{sample}",
        num_inputs = lambda w: len(get_runs_for_sample(w.sample)),
        input_bams = lambda w: " ".join([os.path.join(result_path, "results", sr, "mapped", f"{sr}.filtered.bam") 
                                          for sr in get_runs_for_sample(w.sample)]),
        # Only match if this is a valid sample name (check happens in input function)
        is_valid_sample = lambda w: w.sample in samples,
    resources:
        mem_mb=config.get("mem", "16000"),
    threads: config.get("threads", 2)
    conda:
        "../envs/bowtie2.yaml",
    log:
        "logs/rules/merge_bam_{sample}.log"
    shell:
        """
        mkdir -p $(dirname {output.bam})
        
        if [ {params.num_inputs} -eq 1 ]; then
            # Single run, just copy
            cp {input[0]} {output.bam}
            if [ -f {input[0]}.bai ]; then
                cp {input[0]}.bai {output.bai}
            else
                samtools index {output.bam}
            fi
        else
            # Multiple runs, merge
            samtools merge -f -@ {threads} {output.bam} {params.input_bams} 2> {log}
            samtools index {output.bam} 2>> {log}
        fi
        """
        
# peak calling with MACS2 & samtools and annotation with HOMER
rule peak_calling:
    input:
        bam = os.path.join(result_path,"bam","{sample}", "{sample}.filtered.bam"),
        bai = os.path.join(result_path,"bam","{sample}", "{sample}.filtered.bam.bai"),
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
