# Replicate processing workflow
# For each replicate: BAM -> tagAlign -> pseudo-replicates -> peaks -> combined peaks

# ============================================================================
# Configuration & Constants
# ============================================================================

# Genome references
GENOME_SIZE = config["refs"]["genome_size_bp"]
CHROM_SIZES = config["refs"]["chrom_sizes"]

# Replicate settings
REP_OPTS = config.get("replicates", {})
PVAL_THRESH = REP_OPTS.get("pval_thresh", 1e-3)
SMOOTH_WIN = REP_OPTS.get("smooth_win", 150)
CAP_NUM_PEAK = REP_OPTS.get("cap_num_peak", 300000)
PSEUDO_RANDOM_SEED = REP_OPTS.get("pseudoreplication_random_seed", 0)

# IDR settings
IDR_THRESH = REP_OPTS.get("idr_thresh", 0.05)
IDR_RANK = REP_OPTS.get("idr_rank", "signal.value")

# Calculate shiftsize for MACS2 (negative half of smooth window)
SHIFT_SIZE = -1 * (int(SMOOTH_WIN) // 2)

# Calculate IDR rank column (narrowPeak format)
# Column 7 = signalValue, Column 8 = pValue, Column 9 = qValue
IDR_RANK_COL = {
    "signal.value": "7",
    "p.value": "8",
    "q.value": "9"
}.get(IDR_RANK, "7")

# Calculate negative log10 threshold for IDR filtering
import math
NEG_LOG10_THRESH = f"{-math.log10(float(IDR_THRESH)):.6f}"

# Resource defaults
MEM_MB_DEFAULT = config["resources"].get("mem_mb", 16000)
THREADS_DEFAULT = config["resources"].get("threads", 2)

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
        mem_mb = 4*config["resources"].get("mem_mb", 16000),
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
        bash workflow/scripts/split_pseudoreplicates.sh \
            {input.tagalign} \
            {output.pr1} \
            {output.pr2} \
            {wildcards.replicate} \
            {params.is_paired} \
            {params.random_seed} \
            {log}
        """

# Rule 2: Call peaks on individual sample pseudo-replicates (pr1, pr2)
rule call_peaks_pseudoreplicate:
    input:
        tagalign = os.path.join(result_path, "middle_files", "replicates", "{replicate}", "{replicate}.{pr}.tagAlign.gz"),
    output:
        narrowpeak = os.path.join(result_path, "middle_files", "replicates", "{replicate}", "{replicate}.{pr}.narrowPeak.gz"),
    params:
        out_dir = lambda w: os.path.join(result_path, "middle_files", "replicates", w.replicate),
        prefix = lambda w: os.path.join(result_path, "middle_files", "replicates", w.replicate, f"{w.replicate}.{w.pr}"),
        genome_size = GENOME_SIZE,
        pval_thresh = PVAL_THRESH,
        smooth_win = SMOOTH_WIN,
        cap_num_peak = CAP_NUM_PEAK,
        chrom_sizes = CHROM_SIZES,
        shiftsize = SHIFT_SIZE,
    resources:
        mem_mb = MEM_MB_DEFAULT,
        runtime = 600,
    threads: THREADS_DEFAULT
    wildcard_constraints:
        replicate = "|".join(samples.keys()),
        pr = "pr1|pr2"
    conda:
        "../envs/reproducibility.yaml"
    log:
        "logs/rules/replicates/call_peaks_{replicate}_{pr}.log"
    script:
        "../scripts/call_peaks_macs2.sh"

# Rule 3: Call peaks on full individual sample (for pooled peak reference)
rule call_peaks_individual_sample:
    input:
        tagalign = lambda w: os.path.join(result_path, "middle_files", "bed", f"{w.sample}.tagAlign.gz"),
    output:
        narrowpeak = os.path.join(result_path, "middle_files", "replicates", "{sample}", "{sample}.narrowPeak.gz"),
    params:
        out_dir = lambda w: os.path.join(result_path, "middle_files", "replicates", w.sample),
        prefix = lambda w: os.path.join(result_path, "middle_files", "replicates", w.sample, f"{w.sample}"),
        genome_size = GENOME_SIZE,
        pval_thresh = PVAL_THRESH,
        smooth_win = SMOOTH_WIN,
        cap_num_peak = CAP_NUM_PEAK,
        chrom_sizes = CHROM_SIZES,
        shiftsize = SHIFT_SIZE,
    resources:
        mem_mb = MEM_MB_DEFAULT,
        runtime = 600,
    threads: THREADS_DEFAULT
    wildcard_constraints:
        sample = "|".join(samples.keys())
    conda:
        "../envs/reproducibility.yaml"
    log:
        "logs/rules/replicates/call_peaks_individual_{sample}.log"
    script:
        "../scripts/call_peaks_macs2.sh"

# Rule 4: Pool peaks from individual samples in a group for IDR reference
rule pool_sample_peaks_for_group:
    input:
        peaks = lambda w: [
            os.path.join(result_path, "middle_files", "replicates", sample, f"{sample}.narrowPeak.gz")
            for sample in get_samples_for_replicate(w.group)
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
        idr_peak = os.path.join(result_path, "middle_files", "replicates", "{sample}", "{sample}-pr1_vs_{sample}-pr2.narrowPeak.gz"),
    params:
        out_dir = lambda w: os.path.join(result_path, "middle_files", "replicates", w.sample),
        prefix = lambda w: f"{w.sample}-pr1_vs_{w.sample}-pr2",
        idr_thresh = IDR_THRESH,
        idr_rank = IDR_RANK,
        chrom_sizes = CHROM_SIZES,
        idr_rank_col = IDR_RANK_COL,
        neg_log10_thresh = NEG_LOG10_THRESH,
        method = REP_OPTS.get("peak_combination_method", "idr"),
        nonamecheck = "-nonamecheck" if REP_OPTS.get("nonamecheck", False) else "",
    resources:
        mem_mb = config["resources"].get("mem_mb", 32000),
        runtime = 600,
    threads: THREADS_DEFAULT
    wildcard_constraints:
        sample = "|".join(samples.keys())
    conda:
        "../envs/reproducibility.yaml"
    log:
        "logs/rules/replicates/idr_self_consistency_{sample}.log"
    script:
        "../scripts/idr_analysis.sh"

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
        genome_size = GENOME_SIZE,
        pval_thresh = PVAL_THRESH,
        smooth_win = SMOOTH_WIN,
        cap_num_peak = CAP_NUM_PEAK,
        chrom_sizes = CHROM_SIZES,
        shiftsize = SHIFT_SIZE,
    resources:
        mem_mb = MEM_MB_DEFAULT,
        runtime = 600,
    threads: THREADS_DEFAULT
    conda:
        "../envs/reproducibility.yaml"
    log:
        "logs/rules/replicates/call_peaks_pooled_pr1_{group}.log"
    script:
        "../scripts/call_peaks_macs2.sh"

# Rule 9: Call peaks on pooled-pr2
rule call_peaks_pooled_pr2:
    input:
        tagalign = os.path.join(result_path, "middle_files", "replicates", "groups", "{group}", "{group}.pooled-pr2.tagAlign.gz"),
    output:
        narrowpeak = os.path.join(result_path, "middle_files", "replicates", "groups", "{group}", "{group}.pooled-pr2.narrowPeak.gz"),
    params:
        out_dir = lambda w: os.path.join(result_path, "middle_files", "replicates", "groups", w.group),
        prefix = lambda w: os.path.join(result_path, "middle_files", "replicates", "groups", w.group, f"{w.group}.pooled-pr2"),
        genome_size = GENOME_SIZE,
        pval_thresh = PVAL_THRESH,
        smooth_win = SMOOTH_WIN,
        cap_num_peak = CAP_NUM_PEAK,
        chrom_sizes = CHROM_SIZES,
        shiftsize = SHIFT_SIZE,
    resources:
        mem_mb = MEM_MB_DEFAULT,
        runtime = 600,
    threads: THREADS_DEFAULT
    conda:
        "../envs/reproducibility.yaml"
    log:
        "logs/rules/replicates/call_peaks_pooled_pr2_{group}.log"
    script:
        "../scripts/call_peaks_macs2.sh"

# Rule 10: Call peaks on all pooled tagAligns
rule call_peaks_pooled_all:
    input:
        tagalign = os.path.join(result_path, "middle_files", "replicates", "groups", "{group}", "{group}.pooled-all.tagAlign.gz"),
    output:
        narrowpeak = os.path.join(result_path, "middle_files", "replicates", "groups", "{group}", "{group}.pooled-all.narrowPeak.gz"),
    params:
        out_dir = lambda w: os.path.join(result_path, "middle_files", "replicates", "groups", w.group),
        prefix = lambda w: os.path.join(result_path, "middle_files", "replicates", "groups", w.group, f"{w.group}.pooled-all"),
        genome_size = GENOME_SIZE,
        pval_thresh = PVAL_THRESH,
        smooth_win = SMOOTH_WIN,
        cap_num_peak = CAP_NUM_PEAK,
        chrom_sizes = CHROM_SIZES,
        shiftsize = SHIFT_SIZE,
    resources:
        mem_mb = 4*MEM_MB_DEFAULT,
        runtime = 600,
    threads: THREADS_DEFAULT
    conda:
        "../envs/reproducibility.yaml"
    log:
        "logs/rules/replicates/call_peaks_pooled_all_{group}.log"
    script:
        "../scripts/call_peaks_macs2.sh"

# Rule 11: IDR between pooled-pr1 vs pooled-pr2
rule idr_pooled_pseudoreplicates:
    input:
        pr1_peaks = os.path.join(result_path, "middle_files", "replicates", "groups", "{group}", "{group}.pooled-pr1.narrowPeak.gz"),
        pr2_peaks = os.path.join(result_path, "middle_files", "replicates", "groups", "{group}", "{group}.pooled-pr2.narrowPeak.gz"),
        pooled_peaks = os.path.join(result_path, "middle_files", "replicates", "groups", "{group}", "{group}.pooled-all.narrowPeak.gz"),
    output:
        idr_peak = os.path.join(result_path, "middle_files", "replicates", "groups", "{group}", "{group}.pooled-pr1_vs_pooled-pr2.idr.narrowPeak.gz"),
    params:
        out_dir = lambda w: os.path.join(result_path, "middle_files", "replicates", "groups", w.group),
        prefix = lambda w: f"{w.group}.pooled-pr1_vs_pooled-pr2",
        idr_thresh = IDR_THRESH,
        idr_rank = IDR_RANK,
        chrom_sizes = CHROM_SIZES,
        idr_rank_col = IDR_RANK_COL,
        neg_log10_thresh = NEG_LOG10_THRESH,
        method = REP_OPTS.get("peak_combination_method", "idr"),
        nonamecheck = "-nonamecheck" if REP_OPTS.get("nonamecheck", False) else "",
    resources:
        mem_mb = config["resources"].get("mem_mb", 32000),
        runtime = 600,
    threads: THREADS_DEFAULT
    conda:
        "../envs/reproducibility.yaml"
    log:
        "logs/rules/replicates/idr_pooled_pr_{group}.log"
    script:
        "../scripts/idr_analysis.sh"

# ============================================================================
# PHASE 3: True Replicate Comparisons (pairwise between biological replicates)
# ============================================================================
# Note: Individual sample peaks are created by rule call_peaks_individual_sample (Rule 3)

# Helper function to generate unique replicate pairs (always sorted: rep1 < rep2)
def get_unique_replicate_pairs(group):
    """
    Returns list of unique pairs (rep1, rep2) for a group.
    Enforces rep1 < rep2 lexicographically to ensure unique file paths.
    """
    samples = get_samples_for_replicate(group)
    if len(samples) < 2:
        return []
    
    # Sort samples to ensure deterministic pairing (A vs B, never B vs A)
    samples = sorted(samples)
    
    pairs = []
    for i in range(len(samples)):
        for j in range(i + 1, len(samples)):
            pairs.append((samples[i], samples[j]))
    
    return pairs

# Rule 7: IDR between pairs of true replicates
# Note: This rule assumes rep1 < rep2 lexicographically (enforced upstream via get_unique_replicate_pairs)
rule idr_true_replicates:
    input:
        pr1_peaks = lambda w: os.path.join(result_path, "middle_files", "replicates", w.rep1, f"{w.rep1}.narrowPeak.gz"),
        pr2_peaks = lambda w: os.path.join(result_path, "middle_files", "replicates", w.rep2, f"{w.rep2}.narrowPeak.gz"),
        pooled_peaks = lambda w: os.path.join(result_path, "middle_files", "replicates", "groups", w.group, f"{w.group}.pooled-all.narrowPeak.gz"),
    output:
        idr_peak = os.path.join(result_path, "middle_files", "replicates", "groups", "{group}", "true_replicates", "{rep1}_vs_{rep2}.idr.narrowPeak.gz"),
    params:
        out_dir = lambda w: os.path.join(result_path, "middle_files", "replicates", "groups", w.group, "true_replicates"),
        prefix = lambda w: f"{w.rep1}_vs_{w.rep2}",
        idr_thresh = IDR_THRESH,
        idr_rank = IDR_RANK,
        chrom_sizes = CHROM_SIZES,
        idr_rank_col = IDR_RANK_COL,
        neg_log10_thresh = NEG_LOG10_THRESH,
        method = REP_OPTS.get("peak_combination_method", "idr"),
        nonamecheck = "-nonamecheck" if REP_OPTS.get("nonamecheck", False) else "",
    resources:
        mem_mb = config["resources"].get("mem_mb", 32000),
        runtime = 600,
    threads: THREADS_DEFAULT
    conda:
        "../envs/reproducibility.yaml"
    log:
        "logs/rules/replicates/idr_true_replicates_{group}_{rep1}_vs_{rep2}.log"
    script:
        "../scripts/idr_analysis.sh"

# ============================================================================
# PHASE 4: Replicate Correlation Analysis
# ============================================================================

# Rule: Merge peaks from all replicates in a group for correlation analysis
rule merge_group_peaks_for_correlation:
    input:
        peaks = lambda w: [
            os.path.join(result_path, "middle_files", "replicates", sample, f"{sample}.narrowPeak.gz")
            for sample in get_samples_for_replicate(w.group)
        ],
    output:
        merged_peaks = os.path.join(result_path, "middle_files", "replicates", "groups", "{group}", "{group}.merged_peaks.bed"),
    params:
        out_dir = lambda w: os.path.join(result_path, "middle_files", "replicates", "groups", w.group),
    resources:
        mem_mb = config["resources"].get("mem_mb", 8000),
        runtime = 60,
    threads: 1
    wildcard_constraints:
        group = "|".join(get_replicate_group_ids()) if len(get_replicate_group_ids()) > 0 else "NONE_AVAILABLE"
    conda:
        "../envs/reproducibility.yaml"
    log:
        "logs/rules/replicates/merge_group_peaks_{group}.log"
    shell:
        """
        mkdir -p {params.out_dir}
        # Concatenate all peaks, extract first 4 columns (chr, start, end, name) to create BED format
        # Sort and merge overlapping peaks using bedtools merge
        zcat {input.peaks} | \
        awk 'BEGIN {{OFS="\\t"}} {{print $1, $2, $3, $4}}' | \
        sort -k1,1 -k2,2n | \
        bedtools merge -i - -c 4 -o distinct | \
        awk 'BEGIN {{OFS="\\t"}} {{print $1, $2, $3, "peak_" NR}}' > {output.merged_peaks}
        """

# Rule: Calculate correlation between biological replicates based on read counts in group-merged peaks
rule calculate_replicate_correlation:
    input:
        merged_peaks = lambda w: os.path.join(result_path, "middle_files", "replicates", "groups", w.group, f"{w.group}.merged_peaks.bed"),
        rep1_tagalign = lambda w: os.path.join(result_path, "middle_files", "bed", f"{w.rep1}.tagAlign.gz"),
        rep2_tagalign = lambda w: os.path.join(result_path, "middle_files", "bed", f"{w.rep2}.tagAlign.gz"),
    output:
        tsv = os.path.join(result_path, "middle_files", "replicates", "groups", "{group}", "true_replicates", "{rep1}_vs_{rep2}.correlation.tsv"),
        png = os.path.join(result_path, "middle_files", "replicates", "groups", "{group}", "true_replicates", "{rep1}_vs_{rep2}.correlation.png"),
    params:
        out_dir = lambda w: os.path.join(result_path, "middle_files", "replicates", "groups", w.group, "true_replicates"),
        script = os.path.join(workflow.basedir, "scripts", "calculate_replicate_correlation.py"),
    resources:
        mem_mb = 5 * config["resources"].get("mem_mb", 16000),
        runtime = 300,
    threads: config["resources"].get("threads", 2)
    conda:
        "../envs/reproducibility.yaml"
    wildcard_constraints:
        group = "|".join(get_replicate_group_ids()) if len(get_replicate_group_ids()) > 0 else "NONE_AVAILABLE"
    log:
        "logs/rules/replicates/correlation_{group}_{rep1}_vs_{rep2}.log"
    shell:
        """
        set -euo pipefail

        echo "========================================" >> {log}
        echo "Starting correlation calculation for {wildcards.rep1} vs {wildcards.rep2}" >> {log}
        echo "Group: {wildcards.group}" >> {log}
        echo "========================================" >> {log}

        mkdir -p {params.out_dir}

        # Debug: Check input files exist and have content
        echo "" >> {log}
        echo "--- Step 0: Input file validation ---" >> {log}
        echo "Merged peaks file: {input.merged_peaks}" >> {log}
        if [ -f "{input.merged_peaks}" ]; then
            PEAK_COUNT=$(wc -l < "{input.merged_peaks}" || echo "0")
            echo "  Merged peaks: $PEAK_COUNT lines" >> {log}
            echo "  First 3 lines:" >> {log}
            head -3 "{input.merged_peaks}" >> {log} 2>&1 || echo "  (empty or unreadable)" >> {log}
        else
            echo "  ERROR: Merged peaks file not found!" >> {log}
            exit 1
        fi

        echo "Rep1 tagAlign: {input.rep1_tagalign}" >> {log}
        if [ -f "{input.rep1_tagalign}" ]; then
            REP1_LINES=$(zcat "{input.rep1_tagalign}" | wc -l || echo "0")
            echo "  Rep1 tagAlign: $REP1_LINES lines" >> {log}
        else
            echo "  ERROR: Rep1 tagAlign file not found!" >> {log}
            exit 1
        fi

        echo "Rep2 tagAlign: {input.rep2_tagalign}" >> {log}
        if [ -f "{input.rep2_tagalign}" ]; then
            REP2_LINES=$(zcat "{input.rep2_tagalign}" | wc -l || echo "0")
            echo "  Rep2 tagAlign: $REP2_LINES lines" >> {log}
        else
            echo "  ERROR: Rep2 tagAlign file not found!" >> {log}
            exit 1
        fi

        # 1. Count reads per peak for Rep1 using bedtools coverage
        # CHANGED: Removed -sorted flag and on-the-fly sorting.
        # Bedtools will load the small 'merged_peaks' file into memory and stream the large 'tagAlign'.
        echo "" >> {log}
        echo "--- Step 1: Counting reads for Rep1 ---" >> {log}
        bedtools coverage \
            -a {input.merged_peaks} \
            -b {input.rep1_tagalign} \
            -counts > rep1_counts.tmp 2>> {log}

        if [ ! -s rep1_counts.tmp ]; then
            echo "  ERROR: bedtools coverage produced empty output for Rep1" >> {log}
            exit 1
        fi
        REP1_COUNT_LINES=$(wc -l < rep1_counts.tmp || echo "0")
        REP1_TOTAL=$(awk '{{sum += $5}} END {{print sum+0}}' rep1_counts.tmp || echo "0")
        echo "  Rep1 counts: $REP1_COUNT_LINES peaks, total reads: $REP1_TOTAL" >> {log}
        echo "  First 3 lines of rep1_counts.tmp:" >> {log}
        head -3 rep1_counts.tmp >> {log} 2>&1

        # 2. Count reads per peak for Rep2
        # CHANGED: Removed -sorted flag and on-the-fly sorting.
        echo "" >> {log}
        echo "--- Step 2: Counting reads for Rep2 ---" >> {log}
        bedtools coverage \
            -a {input.merged_peaks} \
            -b {input.rep2_tagalign} \
            -counts > rep2_counts.tmp 2>> {log}

        if [ ! -s rep2_counts.tmp ]; then
            echo "  ERROR: bedtools coverage produced empty output for Rep2" >> {log}
            exit 1
        fi
        REP2_COUNT_LINES=$(wc -l < rep2_counts.tmp || echo "0")
        REP2_TOTAL=$(awk '{{sum += $5}} END {{print sum+0}}' rep2_counts.tmp || echo "0")
        echo "  Rep2 counts: $REP2_COUNT_LINES peaks, total reads: $REP2_TOTAL" >> {log}
        echo "  First 3 lines of rep2_counts.tmp:" >> {log}
        head -3 rep2_counts.tmp >> {log} 2>&1

        # Validate that both count files have the same number of lines
        if [ "$REP1_COUNT_LINES" -ne "$REP2_COUNT_LINES" ]; then
            echo "  WARNING: Rep1 and Rep2 have different number of peaks ($REP1_COUNT_LINES vs $REP2_COUNT_LINES)" >> {log}
        fi

        # 3. Paste results: chr, start, end, name (from rep1), count1, count2
        echo "" >> {log}
        echo "--- Step 3: Combining count files ---" >> {log}
        # rep1_counts.tmp has: chr, start, end, name, count
        # rep2_counts.tmp has: chr, start, end, name, count
        # We want: chr, start, end, name, rep1_count, rep2_count
        paste <(cut -f 1-4 rep1_counts.tmp) <(cut -f 5 rep1_counts.tmp) <(cut -f 5 rep2_counts.tmp) > joint_counts.tmp 2>> {log}

        if [ ! -s joint_counts.tmp ]; then
            echo "  ERROR: paste command produced empty output" >> {log}
            exit 1
        fi
        JOINT_LINES=$(wc -l < joint_counts.tmp || echo "0")
        echo "  Joint counts: $JOINT_LINES lines" >> {log}
        echo "  First 3 lines of joint_counts.tmp:" >> {log}
        head -3 joint_counts.tmp >> {log} 2>&1

        # 4. Run Python script to calculate correlation and generate plot
        echo "" >> {log}
        echo "--- Step 4: Running Python correlation script ---" >> {log}
        echo "  Input: joint_counts.tmp ($JOINT_LINES lines)" >> {log}
        echo "  Output TSV: {output.tsv}" >> {log}
        echo "  Output PNG: {output.png}" >> {log}
        python {params.script} \
            --input joint_counts.tmp \
            --out-tsv {output.tsv} \
            --out-png {output.png} 2>> {log}

        # Verify outputs were created
        echo "" >> {log}
        echo "--- Step 5: Verifying outputs ---" >> {log}
        if [ -f "{output.tsv}" ]; then
            TSV_SIZE=$(wc -l < "{output.tsv}" || echo "0")
            echo "  TSV file created: $TSV_SIZE lines" >> {log}
            cat "{output.tsv}" >> {log} 2>&1 || echo "  (could not read TSV)" >> {log}
        else
            echo "  ERROR: TSV output file not created!" >> {log}
            exit 1
        fi

        if [ -f "{output.png}" ]; then
            PNG_SIZE=$(stat -c%s "{output.png}" 2>/dev/null || stat -f%z "{output.png}" 2>/dev/null || echo "0")
            echo "  PNG file created: $PNG_SIZE bytes" >> {log}
        else
            echo "  ERROR: PNG output file not created!" >> {log}
            exit 1
        fi

        # 5. Clean up temporary files
        echo "" >> {log}
        echo "--- Step 6: Cleaning up temporary files ---" >> {log}
        rm -f rep1_counts.tmp rep2_counts.tmp joint_counts.tmp
        echo "  Temporary files removed" >> {log}

        echo "" >> {log}
        echo "Correlation calculation completed successfully" >> {log}
        echo "========================================" >> {log}
        """


# Helper function to get all correlation TSV files for a group
def get_correlation_pair_files(wildcards):
    """
    Get list of all correlation TSV files for a given group.
    Uses get_unique_replicate_pairs to ensure rep1 < rep2 lexicographically.
    """
    if wildcards.group not in get_replicate_group_ids():
        return []
    
    # Get sorted pairs (rep1 < rep2)
    pairs = get_unique_replicate_pairs(wildcards.group)
    
    pair_files = []
    for rep1, rep2 in pairs:
        pair_file = os.path.join(
            result_path, "middle_files", "replicates", "groups", 
            wildcards.group, "true_replicates", 
            f"{rep1}_vs_{rep2}.correlation.tsv"
        )
        pair_files.append(pair_file)
    
    return pair_files

# Rule: Aggregate correlation results for all pairs in a replicate group
rule aggregate_replicate_correlations:
    input:
        correlation_files = lambda w: get_correlation_pair_files(w),
    output:
        summary_tsv = os.path.join(result_path, "report", "reproducibility", "{group}.replicate_correlations.tsv"),
    params:
        # Pass the specific metadata needed by the script
        pairs = lambda w: [
            (rep1, rep2, os.path.join(result_path, "middle_files", "replicates", "groups", w.group, "true_replicates", f"{rep1}_vs_{rep2}.correlation.tsv"))
            for rep1, rep2 in get_unique_replicate_pairs(w.group)
        ],
    resources:
        mem_mb = 1000,
        runtime = 10,
    threads: 1
    conda:
        "../envs/reproducibility.yaml"
    wildcard_constraints:
        group = "|".join(get_replicate_group_ids()) if len(get_replicate_group_ids()) > 0 else "NONE_AVAILABLE"
    log:
        "logs/rules/replicates/aggregate_correlations_{group}.log"
    script:
        "../scripts/aggregate_correlations.py"

# ============================================================================
# PHASE 5: Reproducibility QC and Optimal Peak Selection
# ============================================================================

# Helper function to get all true replicate pair IDR files for a group
def get_true_replicate_pair_files(wildcards):
    """
    Get list of all true replicate pair IDR files for a given group.
    Uses get_unique_replicate_pairs to ensure rep1 < rep2 lexicographically.
    """
    if wildcards.group not in get_replicate_group_ids():
        return []
    
    # Get sorted pairs (rep1 < rep2)
    pairs = get_unique_replicate_pairs(wildcards.group)
    
    pair_files = []
    for rep1, rep2 in pairs:
        pair_file = os.path.join(
            result_path, "middle_files", "replicates", "groups", 
            wildcards.group, "true_replicates", 
            f"{rep1}_vs_{rep2}.idr.narrowPeak.gz"
        )
        pair_files.append(pair_file)
    
    return pair_files

# Rule 8: Reproducibility QC - Evaluate consistency and select optimal peaks
rule reproducibility_qc:
    input:
        peaks_pr = lambda w: [
            os.path.join(result_path, "middle_files", "replicates", rep, f"{rep}-pr1_vs_{rep}-pr2.narrowPeak.gz")
            for rep in get_samples_for_replicate(w.group)
        ],
        peak_ppr = lambda w: os.path.join(result_path, "middle_files", "replicates", "groups", w.group, f"{w.group}.pooled-pr1_vs_pooled-pr2.idr.narrowPeak.gz"),
        # Ensure all true replicate pairs exist
        true_rep_pairs = lambda w: get_true_replicate_pair_files(w),
    output:
        qc_json = os.path.join(result_path, "report", "reproducibility", "{group}.reproducibility.qc.json"),
        optimal_peak = os.path.join(result_path, "middle_files", "replicates", "groups", "{group}", "{group}.optimal_peak.narrowPeak.gz"),
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


