# Replicate processing workflow
# For each replicate: BAM -> tagAlign -> pseudo-replicates -> peaks -> combined peaks

# Helper function to get replicate sample names from annotation
# Uses replicate_sample_name column to group biological replicates
def get_replicate_names():
    """Get list of replicate group names from sample annotation"""
    if 'replicate_sample_name' not in annot.columns:
        return []
    # Get unique replicate_sample_name values that have >1 sample
    all_replicate_names = annot['replicate_sample_name'].dropna().astype(str).unique().tolist()
    # Only keep replicate names with >1 associated sample
    replicate_names = [rep for rep in all_replicate_names if len(get_reproducibility_sample(rep)) > 1]
    return replicate_names

# Get sample names for a replicate group
def get_samples_for_replicate(replicate_name):
    """Get list of sample names belonging to a replicate group"""
    return get_reproducibility_sample(replicate_name)

# Get BAM files for all samples in a replicate group
def get_replicate_bams(replicate_name):
    """Get BAM file paths for all samples in a replicate group"""
    sample_names = get_samples_for_replicate(replicate_name)
    return [os.path.join(result_path, "important_processed", "bam", f"{sample}.filtered.bam") for sample in sample_names]

# Get tagAlign files for all samples in a replicate group
def get_replicate_tagaligns(replicate_name):
    """Get tagAlign file paths for all samples in a replicate group"""
    sample_names = get_samples_for_replicate(replicate_name)
    return [os.path.join(result_path, "middle_files", "bed", f"{sample}.tagAlign.gz") for sample in sample_names]

# Helper function to get all combined pseudo-replicate peak files
# Returns list of peak files: [rep1-pr1_vs_rep1-pr2.narrowPeak.gz, rep2-pr1_vs_rep2-pr2.narrowPeak.gz, ...]
def get_all_replicate_peaks_pr():
    """Get list of all combined pseudo-replicate peak files for downstream analysis"""
    replicate_names = get_replicate_names()
    return [
        os.path.join(result_path, "important_processed", "replicates", f"{rep}-pr1_vs_{rep}-pr2.narrowPeak.gz")
        for rep in replicate_names
    ]

# Get all replicate group IDs (unique groups, not individual replicates)
def get_replicate_group_ids():
    """Get unique replicate group identifiers"""
    if 'replicate_sample_name' not in annot.columns:
        return []
    # Get unique group names where each group has >1 sample
    group_names = []
    for group in annot['replicate_sample_name'].dropna().astype(str).unique():
        samples_in_group = get_reproducibility_sample(group)
        if len(samples_in_group) > 1:
            group_names.append(group)
    return group_names

# Get all pairs of replicates within a group for pairwise comparisons
def get_replicate_pairs(group_id):
    """Get all pairs (i,j) where i<j within a replicate group"""
    samples_in_group = get_reproducibility_sample(group_id)
    pairs = []
    for i in range(len(samples_in_group)):
        for j in range(i+1, len(samples_in_group)):
            pairs.append((samples_in_group[i], samples_in_group[j]))
    return pairs

# Get all pr1 tagAlign files for a group
def get_group_pr1_tagaligns(group_id):
    """Get all pr1 tagAlign files for replicates in a group"""
    samples_in_group = get_reproducibility_sample(group_id)
    return [
        os.path.join(result_path, "middle_files", "replicates", sample, f"{sample}.pr1.tagAlign.gz")
        for sample in samples_in_group
    ]

# Get all pr2 tagAlign files for a group
def get_group_pr2_tagaligns(group_id):
    """Get all pr2 tagAlign files for replicates in a group"""
    samples_in_group = get_reproducibility_sample(group_id)
    return [
        os.path.join(result_path, "middle_files", "replicates", sample, f"{sample}.pr2.tagAlign.gz")
        for sample in samples_in_group
    ]

# Determine if paired-end based on first sample in replicate group
def is_paired_end(replicate_name):
    """Check if replicate is paired-end (based on first sample)"""
    sample_names = get_samples_for_replicate(replicate_name)
    if not sample_names:
        return False
    # Use first sample's read_type
    first_sample = sample_names[0]
    return samples[first_sample].get("read_type", "single") == "paired"

# Validation function to check if replicate name is valid
def validate_replicate_name(wildcards):
    """Validate that the replicate wildcard corresponds to a valid replicate group"""
    valid_replicates = get_replicate_names()
    if wildcards.replicate not in valid_replicates:
        # Get all sample names for better error message
        all_sample_names = list(samples.keys())
        
        # Check if user mistakenly used a sample_name instead of replicate_sample_name
        if wildcards.replicate in all_sample_names:
            # Find the correct replicate group for this sample
            if 'replicate_sample_name' in annot.columns:
                sample_row = annot[annot['sample_name'] == wildcards.replicate]
                if not sample_row.empty:
                    correct_replicate = sample_row['replicate_sample_name'].iloc[0]
                    raise ValueError(
                        f"ERROR: '{wildcards.replicate}' is a SAMPLE name, not a replicate GROUP name!\n"
                        f"  Sample '{wildcards.replicate}' belongs to replicate group: '{correct_replicate}'\n"
                        f"  Use replicate='{correct_replicate}' instead.\n\n"
                        f"Valid replicate groups (from 'replicate_sample_name' column): {valid_replicates}\n"
                        f"  Each group contains multiple samples for reproducibility analysis."
                    )
        
        raise ValueError(
            f"ERROR: Invalid replicate name '{wildcards.replicate}'!\n"
            f"Valid replicate groups: {valid_replicates}\n"
            f"  These are values from the 'replicate_sample_name' column with >1 sample.\n"
            f"  If you see an empty list, check that:\n"
            f"    1. Your annotation has a 'replicate_sample_name' column\n"
            f"    2. At least one replicate group has >1 sample (needed for reproducibility analysis)"
        )
    return []

# Rule 1: Split each individual sample's tagAlign into pseudo-replicates
rule split_pseudoreplicates_replicate:
    input:
        tagalign = lambda w: os.path.join(result_path, "middle_files", "bed", f"{w.replicate}.tagAlign.gz"),
    output:
        pr1 = os.path.join(result_path, "middle_files", "replicates", "{replicate}", "{replicate}.pr1.tagAlign.gz"),
        pr2 = os.path.join(result_path, "middle_files", "replicates", "{replicate}", "{replicate}.pr2.tagAlign.gz"),
    params:
        out_dir = lambda w: os.path.join(result_path, "middle_files", "replicates", w.replicate),
        is_paired = lambda w: samples[w.replicate].get("read_type", "single") == "paired" if w.replicate in samples else False,
        random_seed = config.get("replicates", {}).get("pseudoreplication_random_seed", 0),
    resources:
        mem_mb = config["resources"].get("mem_mb", 16000),
        runtime = 300,
    threads: 1
    wildcard_constraints:
        replicate = "|".join(samples.keys())
    conda:
        "../envs/pybedtools.yaml"
    log:
        "logs/rules/replicates/split_pseudoreplicates_{replicate}.log"
    shell:
        r"""
        set -euo pipefail

        prefix="{params.out_dir}/{wildcards.replicate}"
        tmp_pr1="$prefix.00"
        tmp_pr2="$prefix.01"

        # Determine random seed
        if [ "{params.random_seed}" = "0" ]; then
            random_seed=$(zcat -f {input.tagalign} | wc -c)
        else
            random_seed="{params.random_seed}"
        fi

        if [ "{params.is_paired}" = "True" ]; then
            # Paired-end: keep pairs together
            nlines=$(zcat -f {input.tagalign} | wc -l)
            nlines=$((nlines / 2))
            nlines=$((nlines / 2 + 1))

            # Each pair = 2 lines. We'll glue them with sed, shuffle, split, and then unglue.
            # Remove ERR trap to prevent SIGPIPE (exit 141) errors from causing hard failure:
            set +e

            zcat -f {input.tagalign} | sed 'N;s/\n/\t/' | \
            shuf --random-source=<(openssl enc -aes-256-ctr -pass pass:$random_seed -nosalt </dev/zero 2>/dev/null) | \
            split -d -l $nlines - "$prefix."

            set -e

            # Restore paired-end format (convert each tabbed line back to 2 lines)
            awk -F $'\t' 'NF==12{{printf "%s\t%s\t%s\t%s\t%s\t%s\n%s\t%s\t%s\t%s\t%s\t%s\n",$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}}' "$tmp_pr1" | gzip -nc > {output.pr1}
            awk -F $'\t' 'NF==12{{printf "%s\t%s\t%s\t%s\t%s\t%s\n%s\t%s\t%s\t%s\t%s\t%s\n",$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}}' "$tmp_pr2" | gzip -nc > {output.pr2}
        else
            # Single-end
            nlines=$(zcat -f {input.tagalign} | wc -l)
            nlines=$((nlines / 2 + 1))

            # Remove ERR trap to prevent SIGPIPE (exit 141) errors from causing hard failure:
            set +e

            zcat -f {input.tagalign} | \
            shuf --random-source=<(openssl enc -aes-256-ctr -pass pass:$random_seed -nosalt </dev/zero 2>/dev/null) | \
            split -d -l $nlines - "$prefix."

            set -e

            gzip -nc "$tmp_pr1" > {output.pr1}
            gzip -nc "$tmp_pr2" > {output.pr2}
        fi

        rm -f "$tmp_pr1" "$tmp_pr2"
        """

# Rule 2: Call peaks on individual sample pseudo-replicates (pr1, pr2) and full sample
rule call_peaks_pseudoreplicate:
    input:
        tagalign = os.path.join(result_path, "middle_files", "replicates", "{replicate}", "{replicate}.{pr}.tagAlign.gz"),
    output:
        narrowpeak = os.path.join(result_path, "middle_files", "replicates", "{replicate}", "{replicate}.{pr}.narrowPeak.gz"),
    params:
        out_dir = lambda w: os.path.join(result_path, "middle_files", "replicates", w.replicate),
        prefix = lambda w: os.path.join(result_path, "middle_files", "replicates", w.replicate, f"{w.replicate}.{w.pr}"),
        genome_size = config["refs"]["genome_size_bp"],
        pval_thresh = config.get("replicates", {}).get("pval_thresh", 1e-3),
        smooth_win = config.get("replicates", {}).get("smooth_win", 150),
        cap_num_peak = config.get("replicates", {}).get("cap_num_peak", 300000),
        chrom_sizes = config["refs"]["chrom_sizes"],
    resources:
        mem_mb = config["resources"].get("mem_mb", 16000),
        runtime = 600,
    threads: config["resources"].get("threads", 2)
    wildcard_constraints:
        replicate = "|".join(samples.keys()),
        pr = "pr1|pr2"
    conda:
        "../envs/reproducibility.yaml"
    log:
        "logs/rules/replicates/call_peaks_{replicate}_{pr}.log"
    shell:
        r"""
        set -euo pipefail

        shiftsize=$(( {params.smooth_win}  / 2 ))
        shiftsize=$((shiftsize * -1))
        
        macs2 callpeak \
            -t {input.tagalign} -f BED -n {params.prefix} -g {params.genome_size} \
            -p {params.pval_thresh} --shift $shiftsize --extsize {params.smooth_win} \
            --nomodel -B --SPMR --keep-dup all --call-summits \
            --outdir {params.out_dir} 2>> {log}

        npeak_tmp="{params.prefix}.tmp"
        npeak_tmp2="{params.prefix}.tmp2"

        LC_COLLATE=C sort -k 8gr,8gr "{params.prefix}_peaks.narrowPeak" | \
        awk 'BEGIN{{OFS="\t"}} {{$4="Peak_"NR}} ($2<0){{$2=0}} ($3<0){{$3=0}} ($10==-1){{$10=int($2+($3-$2+1)/2.0)}} {{print $0}}' > "$npeak_tmp"

        head -n {params.cap_num_peak} "$npeak_tmp" > "$npeak_tmp2"
        bedClip "$npeak_tmp2" {params.chrom_sizes} stdout | gzip -nc > {output.narrowpeak}

        rm -f "$npeak_tmp" "$npeak_tmp2" "{params.prefix}"_*
        """

# Rule 3: Call peaks on full individual sample (for pooled peak reference)
rule call_peaks_individual_sample:
    input:
        tagalign = lambda w: os.path.join(result_path, "middle_files", "bed", f"{w.sample}.tagAlign.gz"),
    output:
        narrowpeak = os.path.join(result_path, "middle_files", "replicates", "{sample}", "{sample}.narrowPeak.gz"),
    params:
        out_dir = lambda w: os.path.join(result_path, "middle_files", "replicates", w.sample),
        prefix = lambda w: os.path.join(result_path, "middle_files", "replicates", w.sample, f"{w.sample}"),
        genome_size = config["refs"]["genome_size_bp"],
        pval_thresh = config.get("replicates", {}).get("pval_thresh", 1e-3),
        smooth_win = config.get("replicates", {}).get("smooth_win", 150),
        cap_num_peak = config.get("replicates", {}).get("cap_num_peak", 300000),
        chrom_sizes = config["refs"]["chrom_sizes"],
    resources:
        mem_mb = config["resources"].get("mem_mb", 16000),
        runtime = 600,
    threads: config["resources"].get("threads", 2)
    wildcard_constraints:
        sample = "|".join(samples.keys())
    conda:
        "../envs/reproducibility.yaml"
    log:
        "logs/rules/replicates/call_peaks_individual_{sample}.log"
    shell:
        r"""
        set -euo pipefail

        mkdir -p {params.out_dir}
        shiftsize=$(( {params.smooth_win}  / 2 ))
        shiftsize=$((shiftsize * -1))
        
        macs2 callpeak \
            -t {input.tagalign} -f BED -n {params.prefix} -g {params.genome_size} \
            -p {params.pval_thresh} --shift $shiftsize --extsize {params.smooth_win} \
            --nomodel -B --SPMR --keep-dup all --call-summits \
            --outdir {params.out_dir} 2>> {log}

        npeak_tmp="{params.prefix}.tmp"
        npeak_tmp2="{params.prefix}.tmp2"

        LC_COLLATE=C sort -k 8gr,8gr "{params.prefix}_peaks.narrowPeak" | \
        awk 'BEGIN{{OFS="\t"}} {{$4="Peak_"NR}} ($2<0){{$2=0}} ($3<0){{$3=0}} ($10==-1){{$10=int($2+($3-$2+1)/2.0)}} {{print $0}}' > "$npeak_tmp"

        head -n {params.cap_num_peak} "$npeak_tmp" > "$npeak_tmp2"
        bedClip "$npeak_tmp2" {params.chrom_sizes} stdout | gzip -nc > {output.narrowpeak}

        rm -f "$npeak_tmp" "$npeak_tmp2" "{params.prefix}"_*
        """

# Rule 4: Pool peaks from individual samples in a group for IDR reference
rule pool_sample_peaks_for_group:
    input:
        peaks = lambda w: [
            os.path.join(result_path, "middle_files", "replicates", sample, f"{sample}.narrowPeak.gz")
            for sample in get_reproducibility_sample(w.group)
        ],
    output:
        pooled_peaks = os.path.join(result_path, "middle_files", "replicates", "groups", "{group}", "{group}.pooled-samples.narrowPeak.gz"),
    params:
        out_dir = lambda w: os.path.join(result_path, "middle_files", "replicates", "groups", w.group),
    resources:
        mem_mb = config["resources"].get("mem_mb", 8000),
        runtime = 60,
    threads: 1
    wildcard_constraints:
        group = "|".join(get_replicate_names()) if len(get_replicate_names()) > 0 else "NONE_AVAILABLE"
    conda:
        "../envs/reproducibility.yaml"
    log:
        "logs/rules/replicates/pool_sample_peaks_{group}.log"
    shell:
        """
        mkdir -p {params.out_dir}
        # Concatenate all peaks, sort, and merge overlapping peaks
        zcat {input.peaks} | sort -k1,1 -k2,2n | gzip -nc > {output.pooled_peaks}
        """

# Rule 5: IDR between pr1 and pr2 for each individual sample (self-consistency)
rule idr_self_consistency:
    input:
        pr1_peaks = os.path.join(result_path, "middle_files", "replicates", "{sample}", "{sample}.pr1.narrowPeak.gz"),
        pr2_peaks = os.path.join(result_path, "middle_files", "replicates", "{sample}", "{sample}.pr2.narrowPeak.gz"),
        pooled_peaks = os.path.join(result_path, "middle_files", "replicates", "{sample}", "{sample}.narrowPeak.gz"),
    output:
        idr_peak = os.path.join(result_path, "important_processed", "replicates", "{sample}-pr1_vs_{sample}-pr2.narrowPeak.gz"),
    params:
        out_dir = lambda w: os.path.join(result_path, "important_processed", "replicates"),
        prefix = lambda w: f"{w.sample}-pr1_vs_{w.sample}-pr2",
        idr_thresh = config.get("replicates", {}).get("idr_thresh", 0.05),
        idr_rank = lambda w: config.get("replicates", {}).get("idr_rank", "signal.value"),
        chrom_sizes = config["refs"]["chrom_sizes"],
        idr_rank_col = lambda w: "8" if config.get("replicates", {}).get("idr_rank", "signal.value") == "signal.value" else "9",
    resources:
        mem_mb = config["resources"].get("mem_mb", 32000),
        runtime = 600,
    threads: config["resources"].get("threads", 2)
    wildcard_constraints:
        sample = "|".join(samples.keys())
    conda:
        "../envs/reproducibility.yaml"
    log:
        "logs/rules/replicates/idr_self_consistency_{sample}.log"
    shell:
        r"""
        mkdir -p {params.out_dir}

        idr_out="{params.out_dir}/{params.prefix}.unthresholded-peaks.txt"
        idr_tmp="{params.out_dir}/{params.prefix}.unthresholded-peaks.txt.tmp"

        idr --samples {input.pr1_peaks} {input.pr2_peaks} --peak-list {input.pooled_peaks} \
            --input-file-type narrowPeak --output-file "$idr_out" \
            --rank {params.idr_rank} --soft-idr-threshold {params.idr_thresh} \
            --plot --use-best-multisummit-IDR --log-output-file {log} 2>&1

        bedClip "$idr_out" {params.chrom_sizes} stdout > "$idr_tmp"

        neg_log10_thresh=$(python3 -c 'import math; print(format(-math.log10(float("{params.idr_thresh}")), ".6f"))')
        rank_col={params.idr_rank_col}

        awk -v thresh="$neg_log10_thresh" 'BEGIN{{OFS="\t"}} $12>=thresh {{
            if ($2<0) $2=0;
            print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12
        }}' "$idr_tmp" | \
        sort -k1,1 -k2,2n | uniq | \
        sort -grk"$rank_col","$rank_col" | \
        awk 'BEGIN{{OFS="\t"}}{{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}}' | \
        gzip -nc > {output.idr_peak}

        rm -f "$idr_out" "$idr_tmp" "{params.out_dir}/{params.prefix}"*.png
        """

# ============================================================================
# PHASE 2: Pooled Pseudo-Replicate Analysis (for replicate groups with >1 sample)
# ============================================================================

# Rule 6: Pool all pr1 tagAligns across all samples in a group
rule pool_group_pr1_tagaligns:
    input:
        pr1_tagaligns = lambda w: get_group_pr1_tagaligns(w.group),
    output:
        pooled_pr1 = os.path.join(result_path, "middle_files", "replicates", "groups", "{group}", "{group}.pooled-pr1.tagAlign.gz"),
    params:
        out_dir = lambda w: os.path.join(result_path, "middle_files", "replicates", "groups", w.group),
    resources:
        mem_mb = config["resources"].get("mem_mb", 16000),
        runtime = 120,
    threads: 1
    conda:
        "../envs/pybedtools.yaml"
    log:
        "logs/rules/replicates/pool_group_pr1_{group}.log"
    shell:
        """
        mkdir -p {params.out_dir}
        zcat {input.pr1_tagaligns} | gzip -nc > {output.pooled_pr1}
        """

# Rule 6: Pool all pr2 tagAligns across all replicates in a group
rule pool_group_pr2_tagaligns:
    input:
        pr2_tagaligns = lambda w: get_group_pr2_tagaligns(w.group),
    output:
        pooled_pr2 = os.path.join(result_path, "middle_files", "replicates", "groups", "{group}", "{group}.pooled-pr2.tagAlign.gz"),
    params:
        out_dir = lambda w: os.path.join(result_path, "middle_files", "replicates", "groups", w.group),
    resources:
        mem_mb = config["resources"].get("mem_mb", 16000),
        runtime = 120,
    threads: 1
    conda:
        "../envs/pybedtools.yaml"
    log:
        "logs/rules/replicates/pool_group_pr2_{group}.log"
    shell:
        """
        mkdir -p {params.out_dir}
        zcat {input.pr2_tagaligns} | gzip -nc > {output.pooled_pr2}
        """

# Rule 7: Pool all replicate tagAligns in a group (for peak calling on all pooled data)
rule pool_group_all_tagaligns:
    input:
        tagaligns = lambda w: get_replicate_tagaligns(w.group),
    output:
        pooled_all = os.path.join(result_path, "middle_files", "replicates", "groups", "{group}", "{group}.pooled-all.tagAlign.gz"),
    params:
        out_dir = lambda w: os.path.join(result_path, "middle_files", "replicates", "groups", w.group),
    resources:
        mem_mb = config["resources"].get("mem_mb", 16000),
        runtime = 120,
    threads: 1
    conda:
        "../envs/pybedtools.yaml"
    log:
        "logs/rules/replicates/pool_group_all_{group}.log"
    shell:
        """
        mkdir -p {params.out_dir}
        zcat {input.tagaligns} | gzip -nc > {output.pooled_all}
        """

# Rule 8: Call peaks on pooled-pr1
rule call_peaks_pooled_pr1:
    input:
        tagalign = os.path.join(result_path, "middle_files", "replicates", "groups", "{group}", "{group}.pooled-pr1.tagAlign.gz"),
    output:
        narrowpeak = os.path.join(result_path, "middle_files", "replicates", "groups", "{group}", "{group}.pooled-pr1.narrowPeak.gz"),
    params:
        out_dir = lambda w: os.path.join(result_path, "middle_files", "replicates", "groups", w.group),
        prefix = lambda w: os.path.join(result_path, "middle_files", "replicates", "groups", w.group, f"{w.group}.pooled-pr1"),
        genome_size = config["refs"]["genome_size_bp"],
        pval_thresh = config.get("replicates", {}).get("pval_thresh", 1e-3),
        smooth_win = config.get("replicates", {}).get("smooth_win", 150),
        cap_num_peak = config.get("replicates", {}).get("cap_num_peak", 300000),
        chrom_sizes = config["refs"]["chrom_sizes"],
    resources:
        mem_mb = config["resources"].get("mem_mb", 16000),
        runtime = 600,
    threads: config["resources"].get("threads", 2)
    conda:
        "../envs/reproducibility.yaml"
    log:
        "logs/rules/replicates/call_peaks_pooled_pr1_{group}.log"
    shell:
        r"""
        set -euo pipefail

        shiftsize=$(( {params.smooth_win}  / 2 ))
        shiftsize=$((shiftsize * -1))
        
        macs2 callpeak \
            -t {input.tagalign} -f BED -n {params.prefix} -g {params.genome_size} \
            -p {params.pval_thresh} --shift $shiftsize --extsize {params.smooth_win} \
            --nomodel -B --SPMR --keep-dup all --call-summits \
            --outdir {params.out_dir} 2>> {log}

        npeak_tmp="{params.prefix}.tmp"
        npeak_tmp2="{params.prefix}.tmp2"

        LC_COLLATE=C sort -k 8gr,8gr "{params.prefix}_peaks.narrowPeak" | \
        awk 'BEGIN{{OFS="\t"}} {{$4="Peak_"NR}} ($2<0){{$2=0}} ($3<0){{$3=0}} ($10==-1){{$10=int($2+($3-$2+1)/2.0)}} {{print $0}}' > "$npeak_tmp"

        head -n {params.cap_num_peak} "$npeak_tmp" > "$npeak_tmp2"
        bedClip "$npeak_tmp2" {params.chrom_sizes} stdout | gzip -nc > {output.narrowpeak}

        rm -f "$npeak_tmp" "$npeak_tmp2" "{params.prefix}"_*
        """

# Rule 9: Call peaks on pooled-pr2
rule call_peaks_pooled_pr2:
    input:
        tagalign = os.path.join(result_path, "middle_files", "replicates", "groups", "{group}", "{group}.pooled-pr2.tagAlign.gz"),
    output:
        narrowpeak = os.path.join(result_path, "middle_files", "replicates", "groups", "{group}", "{group}.pooled-pr2.narrowPeak.gz"),
    params:
        out_dir = lambda w: os.path.join(result_path, "middle_files", "replicates", "groups", w.group),
        prefix = lambda w: os.path.join(result_path, "middle_files", "replicates", "groups", w.group, f"{w.group}.pooled-pr2"),
        genome_size = config["refs"]["genome_size_bp"],
        pval_thresh = config.get("replicates", {}).get("pval_thresh", 1e-3),
        smooth_win = config.get("replicates", {}).get("smooth_win", 150),
        cap_num_peak = config.get("replicates", {}).get("cap_num_peak", 300000),
        chrom_sizes = config["refs"]["chrom_sizes"],
    resources:
        mem_mb = config["resources"].get("mem_mb", 16000),
        runtime = 600,
    threads: config["resources"].get("threads", 2)
    conda:
        "../envs/reproducibility.yaml"
    log:
        "logs/rules/replicates/call_peaks_pooled_pr2_{group}.log"
    shell:
        r"""
        set -euo pipefail

        shiftsize=$(( {params.smooth_win}  / 2 ))
        shiftsize=$((shiftsize * -1))
        
        macs2 callpeak \
            -t {input.tagalign} -f BED -n {params.prefix} -g {params.genome_size} \
            -p {params.pval_thresh} --shift $shiftsize --extsize {params.smooth_win} \
            --nomodel -B --SPMR --keep-dup all --call-summits \
            --outdir {params.out_dir} 2>> {log}

        npeak_tmp="{params.prefix}.tmp"
        npeak_tmp2="{params.prefix}.tmp2"

        LC_COLLATE=C sort -k 8gr,8gr "{params.prefix}_peaks.narrowPeak" | \
        awk 'BEGIN{{OFS="\t"}} {{$4="Peak_"NR}} ($2<0){{$2=0}} ($3<0){{$3=0}} ($10==-1){{$10=int($2+($3-$2+1)/2.0)}} {{print $0}}' > "$npeak_tmp"

        head -n {params.cap_num_peak} "$npeak_tmp" > "$npeak_tmp2"
        bedClip "$npeak_tmp2" {params.chrom_sizes} stdout | gzip -nc > {output.narrowpeak}

        rm -f "$npeak_tmp" "$npeak_tmp2" "{params.prefix}"_*
        """

# Rule 10: Call peaks on all pooled tagAligns
rule call_peaks_pooled_all:
    input:
        tagalign = os.path.join(result_path, "middle_files", "replicates", "groups", "{group}", "{group}.pooled-all.tagAlign.gz"),
    output:
        narrowpeak = os.path.join(result_path, "middle_files", "replicates", "groups", "{group}", "{group}.pooled-all.narrowPeak.gz"),
    params:
        out_dir = lambda w: os.path.join(result_path, "middle_files", "replicates", "groups", w.group),
        prefix = lambda w: os.path.join(result_path, "middle_files", "replicates", "groups", w.group, f"{w.group}.pooled-all"),
        genome_size = config["refs"]["genome_size_bp"],
        pval_thresh = config.get("replicates", {}).get("pval_thresh", 1e-3),
        smooth_win = config.get("replicates", {}).get("smooth_win", 150),
        cap_num_peak = config.get("replicates", {}).get("cap_num_peak", 300000),
        chrom_sizes = config["refs"]["chrom_sizes"],
    resources:
        mem_mb = config["resources"].get("mem_mb", 16000),
        runtime = 600,
    threads: config["resources"].get("threads", 2)
    conda:
        "../envs/reproducibility.yaml"
    log:
        "logs/rules/replicates/call_peaks_pooled_all_{group}.log"
    shell:
        r"""
        set -euo pipefail

        shiftsize=$(( {params.smooth_win}  / 2 ))
        shiftsize=$((shiftsize * -1))
        
        macs2 callpeak \
            -t {input.tagalign} -f BED -n {params.prefix} -g {params.genome_size} \
            -p {params.pval_thresh} --shift $shiftsize --extsize {params.smooth_win} \
            --nomodel -B --SPMR --keep-dup all --call-summits \
            --outdir {params.out_dir} 2>> {log}

        npeak_tmp="{params.prefix}.tmp"
        npeak_tmp2="{params.prefix}.tmp2"

        LC_COLLATE=C sort -k 8gr,8gr "{params.prefix}_peaks.narrowPeak" | \
        awk 'BEGIN{{OFS="\t"}} {{$4="Peak_"NR}} ($2<0){{$2=0}} ($3<0){{$3=0}} ($10==-1){{$10=int($2+($3-$2+1)/2.0)}} {{print $0}}' > "$npeak_tmp"

        head -n {params.cap_num_peak} "$npeak_tmp" > "$npeak_tmp2"
        bedClip "$npeak_tmp2" {params.chrom_sizes} stdout | gzip -nc > {output.narrowpeak}

        rm -f "$npeak_tmp" "$npeak_tmp2" "{params.prefix}"_*
        """

# Rule 11: IDR between pooled-pr1 vs pooled-pr2
rule idr_pooled_pseudoreplicates:
    input:
        pr1_peaks = os.path.join(result_path, "middle_files", "replicates", "groups", "{group}", "{group}.pooled-pr1.narrowPeak.gz"),
        pr2_peaks = os.path.join(result_path, "middle_files", "replicates", "groups", "{group}", "{group}.pooled-pr2.narrowPeak.gz"),
        pooled_peaks = os.path.join(result_path, "middle_files", "replicates", "groups", "{group}", "{group}.pooled-all.narrowPeak.gz"),
    output:
        idr_peak = os.path.join(result_path, "important_processed", "replicates", "groups", "{group}", "{group}.pooled-pr1_vs_pooled-pr2.idr.narrowPeak.gz"),
    params:
        out_dir = lambda w: os.path.join(result_path, "important_processed", "replicates", "groups", w.group),
        prefix = lambda w: f"{w.group}.pooled-pr1_vs_pooled-pr2",
        idr_thresh = config.get("replicates", {}).get("idr_thresh", 0.05),
        idr_rank = lambda w: config.get("replicates", {}).get("idr_rank", "signal.value"),
        chrom_sizes = config["refs"]["chrom_sizes"],
        idr_rank_col = lambda w: "8" if config.get("replicates", {}).get("idr_rank", "signal.value") == "signal.value" else "9",
    resources:
        mem_mb = config["resources"].get("mem_mb", 32000),
        runtime = 600,
    threads: config["resources"].get("threads", 2)
    conda:
        "../envs/reproducibility.yaml"
    log:
        "logs/rules/replicates/idr_pooled_pr_{group}.log"
    shell:
        r"""
        mkdir -p {params.out_dir}

        idr_out="{params.out_dir}/{params.prefix}.unthresholded-peaks.txt"
        idr_tmp="{params.out_dir}/{params.prefix}.unthresholded-peaks.txt.tmp"

        idr --samples {input.pr1_peaks} {input.pr2_peaks} --peak-list {input.pooled_peaks} \
            --input-file-type narrowPeak --output-file "$idr_out" \
            --rank {params.idr_rank} --soft-idr-threshold {params.idr_thresh} \
            --plot --use-best-multisummit-IDR --log-output-file {log} 2>&1

        bedClip "$idr_out" {params.chrom_sizes} stdout > "$idr_tmp"

        neg_log10_thresh=$(python3 -c 'import math; print(format(-math.log10(float("{params.idr_thresh}")), ".6f"))')
        rank_col={params.idr_rank_col}

        awk -v thresh="$neg_log10_thresh" 'BEGIN{{OFS="\t"}} $12>=thresh {{
            if ($2<0) $2=0;
            print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12
        }}' "$idr_tmp" | \
        sort -k1,1 -k2,2n | uniq | \
        sort -grk"$rank_col","$rank_col" | \
        awk 'BEGIN{{OFS="\t"}}{{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}}' | \
        gzip -nc > {output.idr_peak}

        rm -f "$idr_out" "$idr_tmp" "{params.out_dir}/{params.prefix}"*.png
        """

# ============================================================================
# PHASE 3: True Replicate Comparisons (pairwise between biological replicates)
# ============================================================================
# Note: Individual sample peaks are created by rule call_peaks_individual_sample (Rule 3)

# Rule 7: IDR between pairs of true replicates
rule idr_true_replicates:
    input:
        rep1_peaks = lambda w: os.path.join(result_path, "middle_files", "replicates", w.rep1, f"{w.rep1}.narrowPeak.gz"),
        rep2_peaks = lambda w: os.path.join(result_path, "middle_files", "replicates", w.rep2, f"{w.rep2}.narrowPeak.gz"),
        pooled_peaks = lambda w: os.path.join(result_path, "middle_files", "replicates", "groups", w.group, f"{w.group}.pooled-all.narrowPeak.gz"),
    output:
        idr_peak = os.path.join(result_path, "important_processed", "replicates", "groups", "{group}", "true_replicates", "{rep1}_vs_{rep2}.idr.narrowPeak.gz"),
    params:
        out_dir = lambda w: os.path.join(result_path, "important_processed", "replicates", "groups", w.group, "true_replicates"),
        prefix = lambda w: f"{w.rep1}_vs_{w.rep2}",
        idr_thresh = config.get("replicates", {}).get("idr_thresh", 0.05),
        idr_rank = lambda w: config.get("replicates", {}).get("idr_rank", "signal.value"),
        chrom_sizes = config["refs"]["chrom_sizes"],
        idr_rank_col = lambda w: "8" if config.get("replicates", {}).get("idr_rank", "signal.value") == "signal.value" else "9",
    resources:
        mem_mb = config["resources"].get("mem_mb", 32000),
        runtime = 600,
    threads: config["resources"].get("threads", 2)
    conda:
        "../envs/reproducibility.yaml"
    log:
        "logs/rules/replicates/idr_true_replicates_{group}_{rep1}_vs_{rep2}.log"
    shell:
        r"""
        mkdir -p {params.out_dir}

        idr_out="{params.out_dir}/{params.prefix}.unthresholded-peaks.txt"
        idr_tmp="{params.out_dir}/{params.prefix}.unthresholded-peaks.txt.tmp"

        idr --samples {input.rep1_peaks} {input.rep2_peaks} --peak-list {input.pooled_peaks} \
            --input-file-type narrowPeak --output-file "$idr_out" \
            --rank {params.idr_rank} --soft-idr-threshold {params.idr_thresh} \
            --plot --use-best-multisummit-IDR --log-output-file {log} 2>&1

        bedClip "$idr_out" {params.chrom_sizes} stdout > "$idr_tmp"

        neg_log10_thresh=$(python3 -c 'import math; print(format(-math.log10(float("{params.idr_thresh}")), ".6f"))')
        rank_col={params.idr_rank_col}

        awk -v thresh="$neg_log10_thresh" 'BEGIN{{OFS="\t"}} $12>=thresh {{
            if ($2<0) $2=0;
            print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12
        }}' "$idr_tmp" | \
        sort -k1,1 -k2,2n | uniq | \
        sort -grk"$rank_col","$rank_col" | \
        awk 'BEGIN{{OFS="\t"}}{{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}}' | \
        gzip -nc > {output.idr_peak}

        rm -f "$idr_out" "$idr_tmp" "{params.out_dir}/{params.prefix}"*.png
        """

# ============================================================================
# PHASE 4: Reproducibility QC and Optimal Peak Selection
# ============================================================================

# Helper function to get all true replicate pair IDR files for a group
def get_true_replicate_pair_files(wildcards):
    """Get list of all true replicate pair IDR files for a given group"""
    if wildcards.group not in get_replicate_group_ids():
        return []
    
    samples_in_group = get_reproducibility_sample(wildcards.group)
    if len(samples_in_group) < 2:
        return []
    
    pair_files = []
    for i in range(len(samples_in_group)):
        for j in range(i+1, len(samples_in_group)):
            rep_i = samples_in_group[i]
            rep_j = samples_in_group[j]
            pair_file = os.path.join(
                result_path, "important_processed", "replicates", "groups", 
                wildcards.group, "true_replicates", 
                f"{rep_i}_vs_{rep_j}.idr.narrowPeak.gz"
            )
            pair_files.append(pair_file)
    
    return pair_files

# Rule 8: Reproducibility QC - Evaluate consistency and select optimal peaks
rule reproducibility_qc:
    input:
        peaks_pr = lambda w: [
            os.path.join(result_path, "important_processed", "replicates", f"{rep}-pr1_vs_{rep}-pr2.narrowPeak.gz")
            for rep in get_reproducibility_sample(w.group)
        ],
        peak_ppr = lambda w: os.path.join(result_path, "important_processed", "replicates", "groups", w.group, f"{w.group}.pooled-pr1_vs_pooled-pr2.idr.narrowPeak.gz"),
        # Ensure all true replicate pairs exist
        true_rep_pairs = lambda w: get_true_replicate_pair_files(w),
    output:
        qc_json = os.path.join(result_path, "important_processed", "replicates", "groups", "{group}", "{group}.reproducibility.qc.json"),
        optimal_peak = os.path.join(result_path, "important_processed", "replicates", "groups", "{group}", "{group}.optimal_peak.narrowPeak.gz"),
    params:
        prefix = lambda w: w.group,
    resources:
        mem_mb = 4000,
        runtime = 60,
    threads: 1
    conda:
        "../envs/reproducibility.yaml"
    log:
        "logs/rules/replicates/reproducibility_qc_{group}.log"
    script:
        "../scripts/reproducibility_qc.py"


