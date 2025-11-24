prealign_enabled = config["alignment"].get("prealign", {}).get("enabled", True)
prealignments = config["alignment"].get("prealign", {}).get("indices", []) or []
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


def _get_all_trimmed_fastqs_for_sample(sample_name):
    """Get all trimmed FASTQ files for a sample (all runs)"""
    sample_runs = get_runs_for_sample(sample_name)
    r1_files = []
    r2_files = []
    for sr in sample_runs:
        r1_path, r2_path = _trimmed_fastq_paths(sr)
        if r1_path:
            r1_files.append(r1_path)
        if r2_path:
            r2_files.append(r2_path)
    return r1_files, r2_files


def _get_all_prealigned_fastqs_for_sample(sample_name):
    """Get all prealigned FASTQ files for a sample (all runs)"""
    sample_runs = get_runs_for_sample(sample_name)
    r1_files = []
    r2_files = []
    for sr in sample_runs:
        r1_path, r2_path = _prealigned_fastq_paths(sr)
        if r1_path:
            r1_files.append(r1_path)
        if r2_path:
            r2_files.append(r2_path)
    return r1_files, r2_files


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
            aligner = config["alignment"].get("tool", "bowtie2"),
            prealign_dir = lambda w: os.path.join(result_path, "results", w.sample_run, "prealign"),
            is_paired = lambda w: annot.loc[w.sample_run, "read_type"] == "paired",
        resources:
            mem_mb = config["resources"].get("mem_mb", 16000),
        threads: config["resources"].get("threads", 2)
        conda:
            "../envs/bwa.yaml" if config["alignment"].get("tool", "bowtie2") == "bwa-mem2" else "../envs/bowtie2.yaml"
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


if config["alignment"].get("tool", "bowtie2") == "bowtie2":
    # alignment with bowtie2 & samtools (per sample, combining all runs)
    rule align_bowtie2:
        input:
            fasta_fwd = lambda w: (_get_all_prealigned_fastqs_for_sample(w.sample)[0] if has_prealignments else _get_all_trimmed_fastqs_for_sample(w.sample)[0]),
            fasta_rev = lambda w: (_get_all_prealigned_fastqs_for_sample(w.sample)[1] if has_prealignments else _get_all_trimmed_fastqs_for_sample(w.sample)[1]),
            bowtie2_index = os.path.dirname(config["alignment"]["bowtie2"]["index"]),
            adapter_fasta = config["adapters"]["fasta"] if config["adapters"]["fasta"]!="" else [],
            whitelisted_regions = config["refs"]["whitelist"],
        wildcard_constraints:
            sample="|".join(samples.keys())  # Only match actual sample names
        output:
            bam = temp(os.path.join(result_path,"bam","{sample}", "{sample}.filtered.bam")),
            bai = temp(os.path.join(result_path,"bam","{sample}", "{sample}.filtered.bam.bai")),
            bowtie_log = os.path.join(result_path, 'bam', "{sample}", '{sample}.bowtie2.log'),
            bowtie_met = os.path.join(result_path, 'bam', "{sample}", '{sample}.bowtie2.met'),
            samblaster_log = os.path.join(result_path, 'bam', "{sample}", '{sample}.samblaster.log'),
        params:
            sample_name = "{sample}",
            # Format FASTQ files as comma-separated lists for bowtie2
            bowtie2_input = lambda w, input: (
                f"-1 {','.join(input.fasta_fwd)} -2 {','.join(input.fasta_rev)}" 
                if samples[w.sample]["read_type"] == "paired" and input.fasta_rev and len(input.fasta_rev) > 0
                else f"-U {','.join(input.fasta_fwd)}"
            ),
            add_mate_tags = lambda w: "--addMateTags" if samples[w.sample]["read_type"] == "paired" else " ",
            adapter_sequence = "-a " + config["adapters"]["sequence"] if config["adapters"]["sequence"] !="" else " ",
            adapter_fasta = "--adapter_fasta " + config["adapters"]["fasta"] if config["adapters"]["fasta"] !="" else " ",
            sequencing_platform = config["alignment"]["sequencing_platform"],
            bowtie2_index = config["alignment"]["bowtie2"]["index"], # The basename of the index for the reference genome excluding the file endings e.g., *.1.bt2
            bowtie2_local_mode = "--local" if config["alignment"].get("local_mode", False) else "",  # Add this param for local mode
            filtering_flags = lambda w: (
                f"-q 30 -F {config['filtering']['sam_flag']} -f 2 -L {config['refs']['whitelist']}"
                if samples[w.sample]["read_type"] == "paired"
                else f"-q 30 -F {config['filtering']['sam_flag']} -L {config['refs']['whitelist']}"
            ),
        resources:
            mem_mb=config["resources"].get("mem_mb", 16000),
        threads: 4*config["resources"].get("threads", 2)
        conda:
            "../envs/bowtie2.yaml",
        log:
            "logs/rules/align_{sample}.log"
        shell:
            """
            mkdir -p $(dirname {output.bam})
            result_path=$(dirname {output.bam})
            find $result_path -type f -name '*.bam.tmp.*' -exec rm {{}} +;
            
            RG="--rg-id {wildcards.sample} --rg SM:{params.sample_name} --rg PL:{params.sequencing_platform}"

            bowtie2 $RG --very-sensitive --no-discordant -p {threads} --maxins 2000 -x {params.bowtie2_index} \
                {params.bowtie2_local_mode}  \
                --met-file "{output.bowtie_met}" {params.bowtie2_input} 2> "{output.bowtie_log}" | \
                samblaster {params.add_mate_tags} 2> "{output.samblaster_log}" | \
                samtools view {params.filtering_flags} -bhS - 2>> "{output.bowtie_log}" | \
                samtools sort -o "{output.bam}" - 2>> "{output.bowtie_log}";
            
            samtools index "{output.bam}" 2>> "{output.bowtie_log}";
            """

elif config["alignment"].get("tool", "bowtie2") == "bwa-mem2":
    def _sanitize_bwa_min_score_flag():
        value = config["alignment"]["bwa"].get("min_score")
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

    # Using bwa-mem2 only (per sample, combining all runs)

    rule align_bwa_mem:
        input:
            fasta_fwd = lambda w: (_get_all_prealigned_fastqs_for_sample(w.sample)[0] if has_prealignments else _get_all_trimmed_fastqs_for_sample(w.sample)[0]),
            fasta_rev = lambda w: (_get_all_prealigned_fastqs_for_sample(w.sample)[1] if has_prealignments else _get_all_trimmed_fastqs_for_sample(w.sample)[1]),
            # If a custom BWA-MEM2 index is set, use it, otherwise rely on generated genome index files
            index = lambda w: multiext(config["refs"]["fasta"], ".amb", ".ann", ".bwt", ".pac", ".sa", ".0123", ".alt") if not config["alignment"]["bwa"].get("index") or config["alignment"]["bwa"].get("index") in ("", "null", None) else [],
            whitelisted_regions = config["refs"]["whitelist"],
        wildcard_constraints:
            sample="|".join(samples.keys())  # Only match actual sample names
        output:
            bam = temp(os.path.join(result_path,"bam","{sample}", "{sample}.filtered.bam")),
            bai = temp(os.path.join(result_path,"bam","{sample}", "{sample}.filtered.bam.bai")),
            bwa_log = os.path.join(result_path, 'bam', "{sample}", '{sample}.bwa.log'),
            samblaster_log = os.path.join(result_path, 'bam', "{sample}", '{sample}.samblaster.log'),
        params:
            sample_name = "{sample}",
            # Format FASTQ files for bwa-mem2 using process substitution with cat
            # For paired-end: <(cat R1_files...) <(cat R2_files...)
            # For single-end: <(cat R1_files...)
            bwa_input = lambda w, input: (
                f"<(cat {' '.join(input.fasta_fwd)}) <(cat {' '.join(input.fasta_rev)})"
                if samples[w.sample]["read_type"] == "paired" and input.fasta_rev and len(input.fasta_rev) > 0
                else f"<(cat {' '.join(input.fasta_fwd)})"
            ),
            add_mate_tags = lambda w: "--addMateTags" if samples[w.sample]["read_type"] == "paired" else " ",
            sequencing_platform = config["alignment"]["sequencing_platform"],
            bwa_args = config["alignment"]["bwa"].get("extra_args", ""),
            bwa_min_score_flag = _sanitize_bwa_min_score_flag(),
            bwa_m_flag = "-M",  # Mark shorter split hits as secondary
            # Use bwa-mem2 index path if set, otherwise use genome_fasta (for bwa)
            bwa_index_path = lambda w: (config["alignment"]["bwa"].get("index") if config["alignment"]["bwa"].get("index") not in ("", "null", None) else config["refs"]["fasta"]),
            filtering_flags = lambda w: (
                f"-q 30 -F {config['filtering']['sam_flag']} -f 2 -L {config['refs']['whitelist']}"
                if samples[w.sample]["read_type"] == "paired"
                else f"-q 30 -F {config['filtering']['sam_flag']} -L {config['refs']['whitelist']}"
            ),
        resources:
            mem_mb=64000,
            runtime = 800,
        threads: 4*config["resources"].get("threads", 2)
        conda:
            "../envs/bwa.yaml",
        log:
            "logs/rules/align_bwa_mem_{sample}.log"
        shell:
            """
            mkdir -p $(dirname {output.bam})
            result_path=$(dirname {output.bam})
            find $result_path -type f -name '*.bam.tmp.*' -exec rm {{}} +;
            
            RG="@RG\\tID:{wildcards.sample}\\tSM:{params.sample_name}\\tPL:{params.sequencing_platform}"

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
                samtools view {params.filtering_flags} -bhS - 2>> {output.bwa_log} | \
                samtools sort -o {output.bam} - 2>> {output.bwa_log};
            
            samtools index "{output.bam}" 2>> "{output.bwa_log}";
            """

else:
    raise ValueError(f"Unknown aligner: {config['alignment'].get('tool', 'bowtie2')}. Must be 'bowtie2' or 'bwa-mem2'")

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
        # Create BWA index
        bwa-mem2 index {input.fasta} 2> {log}
        """

# Generate alignment stats for per-sample BAM files
rule samtools_process:
    input:
        bam = os.path.join(result_path,"bam","{sample}", "{sample}.filtered.bam"),
        bai = os.path.join(result_path,"bam","{sample}", "{sample}.filtered.bam.bai"),
    output:
        samtools_log = os.path.join(result_path, 'bam', "{sample}", '{sample}.samtools.log'),
        samtools_flagstat_log = os.path.join(result_path, 'bam', "{sample}", '{sample}.samtools_flagstat.log'),
        stats = os.path.join(result_path, 'results', "{sample}", '{sample}.align.stats.tsv'),
    params:
        mitochondria_name = config["refs"].get("mito_name", "chrM"),
    resources:
        mem_mb=config["resources"].get("mem_mb", 16000),
        runtime = 90,
    threads: config["resources"].get("threads", 1)
    conda:
        "../envs/bowtie2.yaml",
    wildcard_constraints:
        sample="|".join(samples.keys())  # Only match actual sample names
    shell:
        '''
            mkdir -p $(dirname {output.stats})
            samtools idxstats "{input.bam}" | awk '{{ sum += $3 + $4; if($1 == "{params.mitochondria_name}") {{ mito_count = $3; }}}}END{{ print "mitochondrial_fraction\t"mito_count/sum }}' > "{output.stats}" 2>> "{output.samtools_log}";
            samtools flagstat "{input.bam}" > "{output.samtools_flagstat_log}" 2>> "{output.samtools_log}";
        '''

# Note: merge_bam rule removed - alignment now directly outputs per-sample BAM files
# combining all runs in a single alignment step
        
# Peak calling with MACS2
rule peak_calling:
    input:
        bam = os.path.join(result_path,"bam","{sample}", "{sample}.filtered.bam"),
        bai = os.path.join(result_path,"bam","{sample}", "{sample}.filtered.bam.bai"),
    output:
        peak_calls = os.path.join(result_path,"results","{sample}","peaks","{sample}_peaks.narrowPeak"),
        macs2_xls = os.path.join(result_path,"results","{sample}","peaks","{sample}_peaks.xls"),
        summits_bed = os.path.join(result_path,"results","{sample}","peaks","{sample}_summits.bed"),
        macs2_log = os.path.join(result_path, 'results', "{sample}", 'peaks', '{sample}.macs2.log'),
    params:
        peaks_dir = os.path.join(result_path,"results","{sample}","peaks"),
        formating = lambda w: '--format BAMPE' if samples["{}".format(w.sample)]["read_type"] == "paired" else '--format BAM',
        genome_size = config["refs"]["genome_size_bp"],
        keep_dup = config["peaks"]["macs2_keep_dup"],
    resources:
        mem_mb=config["resources"].get("mem_mb", 16000),
        runtime = 60,
    threads: config["resources"].get("threads", 2)
    conda:
        "../envs/macs2_homer.yaml",
    log:
        "logs/rules/peak_calling_{sample}.log"
    wildcard_constraints:
        sample="|".join(samples.keys())
    shell:
        """
        mkdir -p {params.peaks_dir}
        
        macs2 callpeak -t {input.bam} {params.formating} \
            --nomodel --keep-dup {params.keep_dup} --extsize 147 -g {params.genome_size} \
            -n {wildcards.sample} \
            --outdir {params.peaks_dir} > "{output.macs2_log}" 2>&1;
        
        # Handle empty peak files - create empty files if MACS2 didn't produce peaks
        if [ ! -f {output.peak_calls} ] || [ ! -s {output.peak_calls} ]; then
            touch {output.peak_calls}
            touch {output.summits_bed}
            echo "Warning: No peaks called for {wildcards.sample}" >> "{output.macs2_log}"
        fi
        """

# Peak annotation and downstream analysis with HOMER
rule peak_annotation:
    input:
        peak_calls = os.path.join(result_path,"results","{sample}","peaks","{sample}_peaks.narrowPeak"),
        summits_bed = os.path.join(result_path,"results","{sample}","peaks","{sample}_summits.bed"),
        bam = os.path.join(result_path,"bam","{sample}", "{sample}.filtered.bam"),
        bai = os.path.join(result_path,"bam","{sample}", "{sample}.filtered.bam.bai"),
        homer_script = os.path.join(HOMER_path,"configureHomer.pl"),
        regulatory_regions = config["refs"]["regulatory_regions"],
    output:
        peak_annot = os.path.join(result_path,"results","{sample}","peaks","{sample}_peaks.narrowPeak.annotated.tsv"),
        peak_annot_log = os.path.join(result_path,"results","{sample}","peaks","{sample}_peaks.narrowPeak.annotated.tsv.log"),
        homer_knownResults = os.path.join(result_path,"results","{sample}","homer","knownResults.txt"),
        homer_log = os.path.join(result_path,"results","{sample}","homer","{sample}.homer.log"),
        stats = os.path.join(result_path, 'results', "{sample}", '{sample}.peak.stats.tsv'),
    params:
        peaks_dir = os.path.join(result_path,"results","{sample}","peaks"),
        homer_dir = os.path.join(result_path,"results","{sample}","homer"),
        homer_bin = os.path.join(HOMER_path,"bin"),
        genome = config["project"]["genome"],
    resources:
        mem_mb=config["resources"].get("mem_mb", 16000),
        runtime = 240,
    threads: config["resources"].get("threads", 2)
    conda:
        "../envs/macs2_homer.yaml",
    log:
        "logs/rules/peak_annotation_{sample}.log"
    wildcard_constraints:
        sample="|".join(samples.keys())
    shell:
        """
        export PATH="{params.homer_bin}:$PATH";
        mkdir -p {params.homer_dir}
        
        # Initialize all output files to ensure they exist
        touch {output.peak_annot}
        touch {output.peak_annot_log}
        touch {output.homer_log}
        touch {output.homer_knownResults}
        touch {output.stats}
        
        # Check if peak file exists and is not empty
        if [ ! -f {input.peak_calls} ] || [ ! -s {input.peak_calls} ]; then
            # Handle empty peak file case
            echo "No peaks found for {wildcards.sample}, creating empty annotation files" > {output.peak_annot_log}
            echo "No peaks found for {wildcards.sample}, skipping HOMER motif finding" > {output.homer_log}
            echo "peaks\t0" > {output.stats}
            echo "frip\t0" >> {output.stats}
            echo "regulatory_fraction\t0" >> {output.stats}
        else
            # Run HOMER annotation
            {params.homer_bin}/annotatePeaks.pl {input.peak_calls} {params.genome} > {output.peak_annot} 2> {output.peak_annot_log};
            
            # Run HOMER motif finding (only if summits file exists and is not empty)
            if [ -f {input.summits_bed} ] && [ -s {input.summits_bed} ]; then
                {params.homer_bin}/findMotifsGenome.pl "{input.summits_bed}" {params.genome} {params.homer_dir} -size 200 -mask > "{output.homer_log}" 2>&1 || echo "HOMER motif finding completed with warnings or errors" >> {output.homer_log}
            else
                echo "No summits file available for motif finding" > {output.homer_log}
            fi
            
            # Calculate statistics
            PEAK_COUNT=$(cat {input.peak_calls} | wc -l)
            echo "peaks\t$PEAK_COUNT" > {output.stats}
            
            TOTAL_READS=$(samtools idxstats {input.bam} | awk '{{sum += $3}}END{{print sum}}');
            
            if [ "$TOTAL_READS" -gt 0 ]; then
                FRIP=$(samtools view -c -L {input.peak_calls} {input.bam} | awk -v total=$TOTAL_READS '{{print $1/total}}')
                echo "frip\t$FRIP" >> {output.stats}
                
                REGULATORY_FRAC=$(samtools view -c -L {input.regulatory_regions} {input.bam} | awk -v total=$TOTAL_READS '{{print $1/total}}')
                echo "regulatory_fraction\t$REGULATORY_FRAC" >> {output.stats}
            else
                echo "frip\t0" >> {output.stats}
                echo "regulatory_fraction\t0" >> {output.stats}
            fi
        fi
        
        # Ensure homer_knownResults exists (findMotifsGenome.pl creates it, but ensure it exists)
        if [ ! -f {output.homer_knownResults} ]; then
            touch {output.homer_knownResults}
        fi
        """

rule merge_peaks:
    input:
        peak_calls = expand(
            os.path.join(result_path, "results", "{sample}", "peaks", "{sample}_peaks.narrowPeak"),
            sample=list(samples.keys()),
        )
    output:
        os.path.join(result_path, "summary", "peaks", "merged_peaks.bed"),
    threads:
        config["resources"].get("threads", 1)
    resources:
        mem_mb=1000,
        runtime = 30,
    params:
        chrom_sizes = config["refs"]["chrom_sizes"],
    log:
        "logs/rules/merge_peaks.log"
    shell:
        """
        mkdir -p $(dirname {output})
        workflow/scripts/merge_peaks --chrom-sizes {params.chrom_sizes} \
            --half-width 250 \
            --output {output} \
            {input.peak_calls}
        """

rule aggregate_stats:
    input:
        peak_stats = os.path.join(result_path, 'results', "{sample}", '{sample}.peak.stats.tsv'),
    output:
        os.path.join(result_path, 'results', "{sample}", '{sample}.stats.tsv'),
    resources:
        mem_mb=config["resources"].get("mem_mb", 1000),
        runtime = 1,
    threads: config["resources"].get("threads", 2)
    log:
        "logs/rules/aggregate_stats_{sample}.log"
    shell:
        """
        cat {input.peak_stats} > {output}
        """


# Generate bigWig tracks from BAM files for visualization
rule tracks:
    input:
        bam = os.path.join(result_path, "bam", "{sample}", "{sample}.filtered.bam"),
        bai = os.path.join(result_path, "bam", "{sample}", "{sample}.filtered.bam.bai"),
    output:
        bw = os.path.join(result_path, "tracks", "{sample}.bw"),
    conda:
        "../envs/dtools.yaml"
    resources:
        mem_mb=config["resources"].get("mem_mb", 16000),
        runtime = 120,
    threads: 8
    wildcard_constraints:
        sample="|".join(samples.keys())  # Only match actual sample names
    log:
        os.path.join("logs","rules","tracks_{sample}.log"),
    shell:
        """
        mkdir -p $(dirname {output.bw})
        bamCoverage -b {input.bam} -o {output.bw} --binSize 10 --smoothLength 50 --normalizeUsing CPM -p {threads} 2> {log}
        """