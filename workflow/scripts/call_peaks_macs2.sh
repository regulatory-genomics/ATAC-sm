#!/bin/bash
# MACS2 peak calling with standard ATAC-seq parameters
# Called by Snakemake - uses snakemake.input, snakemake.output, snakemake.params

set -euo pipefail

# When called via Snakemake script: directive, variables are passed via the environment
TAGALIGN="${snakemake_input[tagalign]}"
OUTPUT_NARROWPEAK="${snakemake_output[narrowpeak]}"
OUT_DIR="${snakemake_params[out_dir]}"
PREFIX="${snakemake_params[prefix]}"
GENOME_SIZE="${snakemake_params[genome_size]}"
PVAL_THRESH="${snakemake_params[pval_thresh]}"
SMOOTH_WIN="${snakemake_params[smooth_win]}"
CAP_NUM_PEAK="${snakemake_params[cap_num_peak]}"
CHROM_SIZES="${snakemake_params[chrom_sizes]}"
SHIFT_SIZE="${snakemake_params[shiftsize]}"

mkdir -p "$OUT_DIR"

# Call peaks with MACS2
macs2 callpeak \
    -t "$TAGALIGN" -f BED -n "$PREFIX" -g "$GENOME_SIZE" \
    -p "$PVAL_THRESH" --shift "$SHIFT_SIZE" --extsize "$SMOOTH_WIN" \
    --nomodel -B --SPMR --keep-dup all --call-summits \
    --outdir "$OUT_DIR"

# Sort peaks by signal value and cap at max number
npeak_tmp="${PREFIX}.tmp"
npeak_tmp2="${PREFIX}.tmp2"

LC_COLLATE=C sort -k 8gr,8gr "${PREFIX}_peaks.narrowPeak" | \
awk 'BEGIN{OFS="\t"} {$4="Peak_"NR} ($2<0){$2=0} ($3<0){$3=0} ($10==-1){$10=int($2+($3-$2+1)/2.0)} {print $0}' > "$npeak_tmp"

head -n "$CAP_NUM_PEAK" "$npeak_tmp" > "$npeak_tmp2"
bedClip "$npeak_tmp2" "$CHROM_SIZES" stdout | gzip -nc > "$OUTPUT_NARROWPEAK"

# Cleanup
rm -f "$npeak_tmp" "$npeak_tmp2" "${PREFIX}"_*

