# Replicate processing workflow
# For each replicate: BAM -> tagAlign -> pseudo-replicates -> peaks -> combined peaks

# Helper function to get replicate names
# Assumes replicates are named as rep1, rep2, etc. or can be configured
def get_replicate_names():
    """Get list of replicate names to process"""
    # Check if replicates are specified in config
    if "replicates" in config and "names" in config["replicates"]:
        return config["replicates"]["names"]
    # Otherwise, try to infer from sample names or annotation
    # For now, return empty list - user should configure replicates
    return []

# Get replicate BAM files
def get_replicate_bam(replicate_name):
    """Get BAM file path for a replicate"""
    # Check if replicate BAM is specified in config
    if "replicates" in config and "bam_files" in config["replicates"]:
        if replicate_name in config["replicates"]["bam_files"]:
            return config["replicates"]["bam_files"][replicate_name]
    # Default: assume replicate BAM follows naming pattern
    return os.path.join(result_path, "important_processed", "bam", f"{replicate_name}.filtered.bam")

# Determine if paired-end based on config or sample annotation
def is_paired_end(replicate_name):
    """Check if replicate is paired-end"""
    if "replicates" in config and "paired_end" in config["replicates"]:
        if isinstance(config["replicates"]["paired_end"], dict):
            return config["replicates"]["paired_end"].get(replicate_name, False)
        return config["replicates"]["paired_end"]
    # Default: check from sample annotation if replicate matches a sample
    for sample_name in samples.keys():
        if replicate_name in sample_name or sample_name in replicate_name:
            return samples[sample_name].get("read_type", "single") == "paired"
    return False

# Helper function to get all combined pseudo-replicate peak files
# Returns list of peak files: [rep1-pr1_vs_rep1-pr2.narrowPeak.gz, rep2-pr1_vs_rep2-pr2.narrowPeak.gz, ...]
def get_all_replicate_peaks_pr():
    """Get list of all combined pseudo-replicate peak files for downstream analysis"""
    replicate_names = get_replicate_names()
    return [
        os.path.join(result_path, "important_processed", "replicates", f"{rep}-pr1_vs_{rep}-pr2.narrowPeak.gz")
        for rep in replicate_names
    ]

# Rule 1: Convert BAM to tagAlign for each replicate
rule bam_to_tagalign_replicate:
    input:
        bam = lambda w: get_replicate_bam(w.replicate),
        bai = lambda w: get_replicate_bam(w.replicate) + ".bai",
    output:
        tagalign = os.path.join(result_path, "middle_files", "replicates", "{replicate}", "{replicate}.tagAlign.gz"),
    params:
        out_dir = lambda w: os.path.join(result_path, "middle_files", "replicates", w.replicate),
        is_paired = lambda w: is_paired_end(w.replicate),
        nth = config["resources"].get("threads", 2),
    resources:
        mem_mb = config["resources"].get("mem_mb", 16000),
        runtime = 300,
    threads: config["resources"].get("threads", 2)
    conda:
        "../envs/pybedtools.yaml"
    log:
        "logs/rules/replicates/bam_to_tagalign_{replicate}.log"
    shell:
        """
        mkdir -p {params.out_dir}
        
        if [ "{params.is_paired}" = "True" ]; then
            # Paired-end: name sort BAM first
            nmsrt_bam="{params.out_dir}/{wildcards.replicate}.nmsrt.bam"
            samtools sort -n -o "$nmsrt_bam" -@ {threads} {input.bam}
            
            # Convert to BEDPE then to tagAlign
            bedpe="{params.out_dir}/{wildcards.replicate}.bedpe.gz"
            bedtools bamtobed -bedpe -mate1 -i "$nmsrt_bam" | gzip -nc > "$bedpe"
            rm -f "$nmsrt_bam"
            
            # Convert BEDPE to tagAlign (two lines per pair)
            zcat -f "$bedpe" | awk 'BEGIN{{OFS="\\t"}}'
                '{{printf "%s\\t%s\\t%s\\tN\\t1000\\t%s\\n%s\\t%s\\t%s\\tN\\t1000\\t%s\\n",'
                '$1,$2,$3,$9,$4,$5,$6,$10}}' | gzip -nc > {output.tagalign}
            rm -f "$bedpe"
        else
            # Single-end: direct conversion
            bedtools bamtobed -i {input.bam} | \
            awk 'BEGIN{{OFS="\\t"}}{{$4="N";$5="1000";print $0}}' | \
            gzip -nc > {output.tagalign}
        fi
        """

# Rule 2: Apply TN5 shift to tagAlign (optional, based on config)
rule tn5_shift_tagalign_replicate:
    input:
        tagalign = os.path.join(result_path, "middle_files", "replicates", "{replicate}", "{replicate}.tagAlign.gz"),
    output:
        shifted_tagalign = os.path.join(result_path, "middle_files", "replicates", "{replicate}", "{replicate}.tn5.tagAlign.gz"),
    params:
        disable_tn5_shift = config.get("replicates", {}).get("disable_tn5_shift", False),
    resources:
        mem_mb = config["resources"].get("mem_mb", 16000),
        runtime = 60,
    threads: 1
    conda:
        "../envs/pybedtools.yaml"
    log:
        "logs/rules/replicates/tn5_shift_{replicate}.log"
    shell:
        """
        if [ "{params.disable_tn5_shift}" = "True" ]; then
            cp {input.tagalign} {output.shifted_tagalign}
        else
            zcat -f {input.tagalign} | awk 'BEGIN {{OFS = "\\t"}} {{
                if ($6 == "+") {{$2 = $2 + 4}} else if ($6 == "-") {{$3 = $3 - 5}} 
                if ($2 >= $3) {{ if ($6 == "+") {{$2 = $3 - 1}} else {{$3 = $2 + 1}} }} 
                print $0}}' | gzip -nc > {output.shifted_tagalign}
        fi
        """

# Rule 3: Split tagAlign into pseudo-replicates
rule split_pseudoreplicates_replicate:
    input:
        tagalign = os.path.join(result_path, "middle_files", "replicates", "{replicate}", "{replicate}.tn5.tagAlign.gz"),
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
        """
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
            
            zcat -f {input.tagalign} | sed 'N;s/\\n/\\t/' | \
            shuf --random-source=<(openssl enc -aes-256-ctr -pass pass:$random_seed -nosalt </dev/zero 2>/dev/null) | \
            split -d -l $nlines - "$prefix."
            
            # Restore paired-end format
            zcat -f "$tmp_pr1" | awk 'BEGIN{{OFS="\\t"}} '
                '{{printf "%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n",'
                '$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}}' | gzip -nc > {output.pr1}
            
            zcat -f "$tmp_pr2" | awk 'BEGIN{{OFS="\\t"}} '
                '{{printf "%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n",'
                '$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}}' | gzip -nc > {output.pr2}
        else
            # Single-end
            nlines=$(zcat -f {input.tagalign} | wc -l)
            nlines=$((nlines / 2 + 1))
            
            zcat -f {input.tagalign} | \
            shuf --random-source=<(openssl enc -aes-256-ctr -pass pass:$random_seed -nosalt </dev/zero 2>/dev/null) | \
            split -d -l $nlines - "$prefix."
            
            gzip -nc "$tmp_pr1" > {output.pr1}
            gzip -nc "$tmp_pr2" > {output.pr2}
        fi
        
        rm -f "$tmp_pr1" "$tmp_pr2"
        """

# Rule 4: Call peaks on pseudo-replicate 1
rule call_peaks_pr1_replicate:
    input:
        tagalign = os.path.join(result_path, "middle_files", "replicates", "{replicate}", "{replicate}.pr1.tagAlign.gz"),
    output:
        narrowpeak = os.path.join(result_path, "middle_files", "replicates", "{replicate}", "{replicate}.pr1.narrowPeak.gz"),
    params:
        out_dir = lambda w: os.path.join(result_path, "middle_files", "replicates", w.replicate),
        prefix = lambda w: os.path.join(result_path, "middle_files", "replicates", w.replicate, w.replicate + ".pr1"),
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
        "../envs/macs2_homer.yaml"
    log:
        "logs/rules/replicates/call_peaks_pr1_{replicate}.log"
    shell:
        """
        shiftsize=$(echo "{params.smooth_win} / 2" | bc | awk '{{printf "%.0f", $1}}')
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
        awk 'BEGIN{{OFS="\\t"}}{{'
            '$4="Peak_"NR; if ($2<0) $2=0; if ($3<0) $3=0; if ($10==-1) '
            '$10=$2+int(($3-$2+1)/2.0); print $0}}' > "$npeak_tmp"
        
        head -n {params.cap_num_peak} "$npeak_tmp" > "$npeak_tmp2"
        
        # Clip peaks between 0-chromSize
        bedClip "$npeak_tmp2" {params.chrom_sizes} stdout | gzip -nc > {output.narrowpeak}
        
        rm -f "$npeak_tmp" "$npeak_tmp2" "{params.prefix}"_*
        """

# Rule 5: Call peaks on pseudo-replicate 2
rule call_peaks_pr2_replicate:
    input:
        tagalign = os.path.join(result_path, "middle_files", "replicates", "{replicate}", "{replicate}.pr2.tagAlign.gz"),
    output:
        narrowpeak = os.path.join(result_path, "middle_files", "replicates", "{replicate}", "{replicate}.pr2.narrowPeak.gz"),
    params:
        out_dir = lambda w: os.path.join(result_path, "middle_files", "replicates", w.replicate),
        prefix = lambda w: os.path.join(result_path, "middle_files", "replicates", w.replicate, w.replicate + ".pr2"),
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
        "../envs/macs2_homer.yaml"
    log:
        "logs/rules/replicates/call_peaks_pr2_{replicate}.log"
    shell:
        """
        shiftsize=$(echo "{params.smooth_win} / 2" | bc | awk '{{printf "%.0f", $1}}')
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
        awk 'BEGIN{{OFS="\\t"}}{{'
            '$4="Peak_"NR; if ($2<0) $2=0; if ($3<0) $3=0; if ($10==-1) '
            '$10=$2+int(($3-$2+1)/2.0); print $0}}' > "$npeak_tmp"
        
        head -n {params.cap_num_peak} "$npeak_tmp" > "$npeak_tmp2"
        
        # Clip peaks between 0-chromSize
        bedClip "$npeak_tmp2" {params.chrom_sizes} stdout | gzip -nc > {output.narrowpeak}
        
        rm -f "$npeak_tmp" "$npeak_tmp2" "{params.prefix}"_*
        """

# Rule 6: Combine pseudo-replicate peaks using IDR or overlap
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
        # Convert IDR rank to column number
        idr_rank_col = lambda w: "8" if config.get("replicates", {}).get("idr_rank", "signal.value") == "signal.value" else "9",
    resources:
        mem_mb = config["resources"].get("mem_mb", 32000),
        runtime = 600,
    threads: config["resources"].get("threads", 2)
    conda:
        "../envs/macs2_homer.yaml"
    log:
        "logs/rules/replicates/combine_peaks_{replicate}.log"
    shell:
        """
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
            
            # Filter by IDR threshold and convert to narrowPeak format
            neg_log10_thresh=$(echo "-l({params.idr_thresh})/l(10)" | bc -l | awk '{{printf "%.6f", $1}}')
            rank_col={params.idr_rank_col}
            
            awk -v thresh="$neg_log10_thresh" 'BEGIN{{OFS="\\t"}} $12>=thresh {{
                if ($2<0) $2=0; 
                print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}}' "$idr_tmp" | \
            sort -k1,1 -k2,2n | uniq | \
            sort -grk$rank_col,$rank_col | \
            awk 'BEGIN{{OFS="\\t"}}{{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}}' | \
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
            awk 'BEGIN{{FS="\\t";OFS="\\t"}} {{
                s1=$3-$2; s2=$13-$12; 
                if (($21/s1 >= 0.5) || ($21/s2 >= 0.5)) {{print $0}}
            }}' | \
            cut -f 1-10 | sort -k1,1 -k2,2n | uniq | \
            # Second intersection: result vs pr2
            intersectBed {params.nonamecheck} -wo -a stdin -b "$tmp2" | \
            awk 'BEGIN{{FS="\\t";OFS="\\t"}} {{
                s1=$3-$2; s2=$13-$12; 
                if (($21/s1 >= 0.5) || ($21/s2 >= 0.5)) {{print $0}}
            }}' | \
            cut -f 1-10 | sort -k1,1 -k2,2n | uniq | gzip -nc > {output.combined_peaks}
            
            rm -f "$tmp1" "$tmp2" "$tmp_pooled"
        fi
        """

# Rule 7: Create pooled peaks for IDR (peak calling on full replicate)
rule call_pooled_peaks_replicate:
    input:
        tagalign = os.path.join(result_path, "middle_files", "replicates", "{replicate}", "{replicate}.tn5.tagAlign.gz"),
    output:
        pooled_peaks = os.path.join(result_path, "middle_files", "replicates", "{replicate}", "{replicate}.pooled.narrowPeak.gz"),
    params:
        out_dir = lambda w: os.path.join(result_path, "middle_files", "replicates", w.replicate),
        prefix = lambda w: os.path.join(result_path, "middle_files", "replicates", w.replicate, w.replicate + ".pooled"),
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
        "../envs/macs2_homer.yaml"
    log:
        "logs/rules/replicates/call_pooled_peaks_{replicate}.log"
    shell:
        """
        shiftsize=$(echo "{params.smooth_win} / 2" | bc | awk '{{printf "%.0f", $1}}')
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
        awk 'BEGIN{{OFS="\\t"}}{{'
            '$4="Peak_"NR; if ($2<0) $2=0; if ($3<0) $3=0; if ($10==-1) '
            '$10=$2+int(($3-$2+1)/2.0); print $0}}' > "$npeak_tmp"
        
        head -n {params.cap_num_peak} "$npeak_tmp" > "$npeak_tmp2"
        
        # Clip peaks between 0-chromSize
        bedClip "$npeak_tmp2" {params.chrom_sizes} stdout | gzip -nc > {output.pooled_peaks}
        
        rm -f "$npeak_tmp" "$npeak_tmp2" "{params.prefix}"_*
        """

