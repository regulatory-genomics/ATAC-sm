#!/bin/bash
# Convert BAM to tagAlign format with TN5 shift
# Usage: bam_to_tagalign.sh <bam> <output> <sample> <is_paired> <disable_tn5_shift> <threads> <bed_dir> <log>

set -euo pipefail

# Parse arguments
BAM="$1"
OUTPUT="$2"
SAMPLE="$3"
IS_PAIRED="$4"
DISABLE_TN5_SHIFT="$5"
THREADS="$6"
BED_DIR="$7"
LOG="$8"

# Create directories
mkdir -p "$BED_DIR"
TMP_DIR="$BED_DIR/tmp_$SAMPLE"
mkdir -p "$TMP_DIR"

# Cleanup function
cleanup() {
    rm -rf "$TMP_DIR"
}
trap cleanup EXIT

TAGALIGN_TMP="$TMP_DIR/$SAMPLE.tagAlign.gz"

if [ "$IS_PAIRED" = "True" ]; then
    # Paired-end: name sort BAM
    NMSRT_BAM="$TMP_DIR/$SAMPLE.nmsrt.bam"
    samtools sort -n -o "$NMSRT_BAM" -@ "$THREADS" "$BAM" 2>> "$LOG"
    
    # BEDPE: bamtobed -bedpe, then to tagAlign (both mates)
    BEDPE="$TMP_DIR/$SAMPLE.bedpe.gz"
    bedtools bamtobed -bedpe -mate1 -i "$NMSRT_BAM" 2>> "$LOG" | gzip -nc > "$BEDPE"
    rm -f "$NMSRT_BAM"
    
    # Convert BEDPE to tagAlign (each mate as a separate tag)
    zcat -f "$BEDPE" | awk 'BEGIN{OFS="\t"} { if($1!="#" && NF>=10) {printf "%s\t%s\t%s\tN\t1000\t%s\n%s\t%s\t%s\tN\t1000\t%s\n", $1,$2,$3,$9,$4,$5,$6,$10} }' | gzip -nc > "$TAGALIGN_TMP"
    rm -f "$BEDPE"
else
    # Single-end: direct conversion
    bedtools bamtobed -i "$BAM" 2>> "$LOG" | awk 'BEGIN{OFS="\t"}{ $4="N"; $5="1000"; print $0 }' | gzip -nc > "$TAGALIGN_TMP"
fi

# TN5 shift
if [ "$DISABLE_TN5_SHIFT" = "True" ]; then
    cp "$TAGALIGN_TMP" "$OUTPUT"
else
    zcat -f "$TAGALIGN_TMP" | awk 'BEGIN{OFS="\t"} { if($6=="+"){$2=$2+4} else if($6=="-"){$3=$3-5} if($2>=$3){if($6=="+"){$2=$3-1} else{$3=$2+1}} print $0 }' | gzip -nc > "$OUTPUT"
fi

# Cleanup handled by trap EXIT

