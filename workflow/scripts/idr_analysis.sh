#!/bin/bash
# IDR or naive overlap analysis between two peak files
# Called by Snakemake - uses snakemake.input, snakemake.output, snakemake.params

set -euo pipefail

# When called via Snakemake script: directive, variables are passed via the environment
PR1_PEAKS="${snakemake_input[pr1_peaks]}"
PR2_PEAKS="${snakemake_input[pr2_peaks]}"
POOLED_PEAKS="${snakemake_input[pooled_peaks]}"
OUTPUT_PEAK="${snakemake_output[idr_peak]}"
OUT_DIR="${snakemake_params[out_dir]}"
PREFIX="${snakemake_params[prefix]}"
IDR_THRESH="${snakemake_params[idr_thresh]}"
IDR_RANK="${snakemake_params[idr_rank]}"
CHROM_SIZES="${snakemake_params[chrom_sizes]}"
RANK_COL="${snakemake_params[idr_rank_col]}"
NEG_LOG10_THRESH="${snakemake_params[neg_log10_thresh]}"
METHOD="${snakemake_params[method]}"
NONAMECHECK="${snakemake_params[nonamecheck]}"

mkdir -p "$OUT_DIR"

if [ "$METHOD" = "idr" ]; then
    # IDR method
    IDR_OUT="${OUT_DIR}/${PREFIX}.idr${IDR_THRESH}.unthresholded-peaks.txt"
    IDR_TMP="${OUT_DIR}/${PREFIX}.idr${IDR_THRESH}.unthresholded-peaks.txt.tmp"

    # Run IDR
    idr --samples "$PR1_PEAKS" "$PR2_PEAKS" --peak-list "$POOLED_PEAKS" \
        --input-file-type narrowPeak --output-file "$IDR_OUT" \
        --rank "$IDR_RANK" --soft-idr-threshold "$IDR_THRESH" \
        --plot --use-best-multisummit-IDR

    # Clip peaks to chromosome sizes
    bedClip "$IDR_OUT" "$CHROM_SIZES" stdout > "$IDR_TMP"

    # Filter by IDR threshold and convert to narrowPeak format
    awk -v thresh="$NEG_LOG10_THRESH" 'BEGIN{OFS="\t"} $12>=thresh {
        if ($2<0) $2=0;
        print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12
    }' "$IDR_TMP" | \
    sort -k1,1 -k2,2n | uniq | \
    sort -grk"${RANK_COL}","${RANK_COL}" | \
    awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}' | \
    gzip -nc > "$OUTPUT_PEAK"

    # Cleanup
    rm -f "$IDR_OUT" "$IDR_TMP" "${OUT_DIR}/${PREFIX}"*.png
else
    # Naive overlap method
    TMP1="${OUT_DIR}/${PREFIX}.tmp1.narrowPeak"
    TMP2="${OUT_DIR}/${PREFIX}.tmp2.narrowPeak"
    TMP_POOLED="${OUT_DIR}/${PREFIX}.tmp_pooled.narrowPeak"

    zcat -f "$PR1_PEAKS" > "$TMP1"
    zcat -f "$PR2_PEAKS" > "$TMP2"
    zcat -f "$POOLED_PEAKS" > "$TMP_POOLED"

    intersectBed $NONAMECHECK -wo -a "$TMP_POOLED" -b "$TMP1" | \
    awk 'BEGIN{FS="\t";OFS="\t"} {
        s1=$3-$2; s2=$13-$12;
        if (($21/s1 >= 0.5) || ($21/s2 >= 0.5)) {print $0}
    }' | \
    cut -f 1-10 | sort -k1,1 -k2,2n | uniq | \
    intersectBed $NONAMECHECK -wo -a stdin -b "$TMP2" | \
    awk 'BEGIN{FS="\t";OFS="\t"} {
        s1=$3-$2; s2=$13-$12;
        if (($21/s1 >= 0.5) || ($21/s2 >= 0.5)) {print $0}
    }' | \
    cut -f 1-10 | sort -k1,1 -k2,2n | uniq | gzip -nc > "$OUTPUT_PEAK"

    rm -f "$TMP1" "$TMP2" "$TMP_POOLED"
fi

