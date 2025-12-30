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
# PHASE 4: Reproducibility QC and Optimal Peak Selection
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


