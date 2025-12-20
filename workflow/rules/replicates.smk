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

# Determine if paired-end based on first sample in replicate group
def is_paired_end(replicate_name):
    """Check if replicate is paired-end (based on first sample)"""
    sample_names = get_samples_for_replicate(replicate_name)
    if not sample_names:
        return False
    # Use first sample's read_type
    first_sample = sample_names[0]
    return samples[first_sample].get("read_type", "single") == "paired"

# Rule 1: Pool tagAlign files from all samples in a replicate group
rule pool_replicate_tagaligns:
    input:
        tagaligns = lambda w: get_replicate_tagaligns(w.replicate),
    output:
        pooled_tagalign = os.path.join(result_path, "middle_files", "replicates", "{replicate}", "{replicate}.pooled.tagAlign.gz"),
    params:
        out_dir = lambda w: os.path.join(result_path, "middle_files", "replicates", w.replicate),
    resources:
        mem_mb = config["resources"].get("mem_mb", 16000),
        runtime = 60,
    threads: 1
    conda:
        "../envs/pybedtools.yaml"
    log:
        "logs/rules/replicates/pool_tagaligns_{replicate}.log"
    shell:
        """
        mkdir -p {params.out_dir}
        # Pool all tagAlign files by concatenating
        zcat {input.tagaligns} | gzip -nc > {output.pooled_tagalign}
        """


# Rule 2: Split pooled tagAlign into pseudo-replicates
rule split_pseudoreplicates_replicate:
    input:
        tagalign = os.path.join(result_path, "middle_files", "replicates", "{replicate}", "{replicate}.pooled.tagAlign.gz"),
    output:
        pr1 = os.path.join(result_path, "middle_files", "replicates", "{replicate}", "{replicate}.pr1.tagAlign.gz"),
        pr2 = os.path.join(result_path, "middle_files", "replicates", "{replicate}", "{replicate}.pr2.tagAlign.gz"),
    params:
        out_dir = lambda w: os.path.join(result_path, "middle_files", "replicates", w.replicate),
        is_paired = lambda w: is_paired_end(w.replicate),
        random_seed = config.get("replicates", {}).get("pseudoreplication_random_seed", 0),
    resources:
        mem_mb = config["resources"].get("mem_mb", 16000),
        runtime = 300,
    threads: 1
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

# Rule 3: Call peaks on pseudo-replicates (handles both pr1 and pr2)
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
        pr="pr1|pr2|pooled"
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

        # Sort and cap peaks
        npeak_tmp="{params.prefix}.tmp"
        npeak_tmp2="{params.prefix}.tmp2"

        LC_COLLATE=C sort -k 8gr,8gr "{params.prefix}_peaks.narrowPeak" | \
        awk 'BEGIN{{OFS="\t"}} {{$4="Peak_"NR}} ($2<0){{$2=0}} ($3<0){{$3=0}} ($10==-1){{$10=int($2+($3-$2+1)/2.0)}} {{print $0}}' > "$npeak_tmp"

        head -n {params.cap_num_peak} "$npeak_tmp" > "$npeak_tmp2"

        # Clip peaks between 0-chromSize
        bedClip "$npeak_tmp2" {params.chrom_sizes} stdout | gzip -nc > {output.narrowpeak}

        rm -f "$npeak_tmp" "$npeak_tmp2" "{params.prefix}"_*
        """

# Rule 4: Combine pseudo-replicate peaks using IDR or naive overlap
rule combine_pseudoreplicate_peaks:
    input:
        pr1_peaks = os.path.join(result_path, "middle_files", "replicates", "{replicate}", "{replicate}.pr1.narrowPeak.gz"),
        pr2_peaks = os.path.join(result_path, "middle_files", "replicates", "{replicate}", "{replicate}.pr2.narrowPeak.gz"),
        pooled_peaks = os.path.join(result_path, "middle_files", "replicates", "{replicate}", "{replicate}.pooled.narrowPeak.gz"),
    output:
        combined_peaks = os.path.join(result_path, "important_processed", "replicates", "{replicate}-pr1_vs_{replicate}-pr2.narrowPeak.gz"),
    params:
        out_dir = lambda w: os.path.join(result_path, "important_processed", "replicates"),
        method = lambda w: config.get("replicates", {}).get("peak_combination_method", "idr"),  # "idr" or "overlap"
        prefix = lambda w: w.replicate + "-pr1_vs_" + w.replicate + "-pr2",
        # IDR parameters
        idr_thresh = config.get("replicates", {}).get("idr_thresh", 0.05),
        idr_rank = lambda w: config.get("replicates", {}).get("idr_rank", "signal.value"),  # "signal.value" (col 8) or "p.value" (col 9)
        # Overlap parameters
        nonamecheck = lambda w: "-nonamecheck" if config.get("replicates", {}).get("nonamecheck", False) else "",
        chrom_sizes = config["refs"]["chrom_sizes"],
        blacklist = config.get("refs", {}).get("blacklist", ""),
        idr_rank_col = lambda w: "8" if config.get("replicates", {}).get("idr_rank", "signal.value") == "signal.value" else "9",
    resources:
        mem_mb = config["resources"].get("mem_mb", 32000),
        runtime = 600,
    threads: config["resources"].get("threads", 2)
    conda:
        "../envs/reproducibility.yaml"
    log:
        "logs/rules/replicates/combine_peaks_{replicate}.log"
    shell:
        r"""
        mkdir -p {params.out_dir}
        method="{params.method}"

        if [ "$method" = "idr" ]; then
            # Use IDR to combine peaks
            idr_out="{params.out_dir}/{params.prefix}.idr{params.idr_thresh}.unthresholded-peaks.txt"
            idr_tmp="{params.out_dir}/{params.prefix}.idr{params.idr_thresh}.unthresholded-peaks.txt.tmp"
            idr_12col="{params.out_dir}/{params.prefix}.idr{params.idr_thresh}.narrowPeak.12-col.bed.gz"

            idr --samples {input.pr1_peaks} {input.pr2_peaks} --peak-list {input.pooled_peaks} \
                --input-file-type narrowPeak --output-file "$idr_out" \
                --rank {params.idr_rank} --soft-idr-threshold {params.idr_thresh} \
                --plot --use-best-multisummit-IDR --log-output-file {log} 2>&1

            # Clip peaks
            bedClip "$idr_out" {params.chrom_sizes} stdout > "$idr_tmp"

            # Filter by IDR threshold and convert to narrowPeak format.
            # Use column 12 ("$12") for IDR score (neglog10 IDR score)
            neg_log10_thresh=$(python3 -c 'import math; print(format(-math.log10(float("{params.idr_thresh}")), ".6f"))')

            rank_col={params.idr_rank_col}

            awk -v thresh="$neg_log10_thresh" 'BEGIN{{OFS="\t"}} $12>=thresh {{
                if ($2<0) $2=0;
                print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12
            }}' "$idr_tmp" | \
            sort -k1,1 -k2,2n | uniq | \
            sort -grk"$rank_col","$rank_col" | \
            awk 'BEGIN{{OFS="\t"}}{{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}}' | \
            gzip -nc > {output.combined_peaks}

            rm -f "$idr_out" "$idr_tmp" "$idr_12col" "{params.out_dir}/{params.prefix}"*.png
        else
            # Use naive overlap to combine peaks
            tmp1="{params.out_dir}/{params.prefix}.tmp1.narrowPeak"
            tmp2="{params.out_dir}/{params.prefix}.tmp2.narrowPeak"
            tmp_pooled="{params.out_dir}/{params.prefix}.tmp_pooled.narrowPeak"

            zcat -f {input.pr1_peaks} > "$tmp1"
            zcat -f {input.pr2_peaks} > "$tmp2"
            zcat -f {input.pooled_peaks} > "$tmp_pooled"

            # Find overlapping peaks (>=50% overlap)
            # First intersection: pooled vs pr1
            intersectBed {params.nonamecheck} -wo -a "$tmp_pooled" -b "$tmp1" | \
            awk 'BEGIN{{FS="\t";OFS="\t"}} {{
                s1=$3-$2; s2=$13-$12;
                if (($21/s1 >= 0.5) || ($21/s2 >= 0.5)) {{print $0}}
            }}' | \
            cut -f 1-10 | sort -k1,1 -k2,2n | uniq | \
            # Second intersection: result vs pr2
            intersectBed {params.nonamecheck} -wo -a stdin -b "$tmp2" | \
            awk 'BEGIN{{FS="\t";OFS="\t"}} {{
                s1=$3-$2; s2=$13-$12;
                if (($21/s1 >= 0.5) || ($21/s2 >= 0.5)) {{print $0}}
            }}' | \
            cut -f 1-10 | sort -k1,1 -k2,2n | uniq | gzip -nc > {output.combined_peaks}

            rm -f "$tmp1" "$tmp2" "$tmp_pooled"
        fi
        """
