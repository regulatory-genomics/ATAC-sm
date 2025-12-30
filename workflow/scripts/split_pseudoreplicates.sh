#!/bin/bash
# Split tagAlign file into two pseudo-replicates (pr1 and pr2)
# Usage: split_pseudoreplicates.sh <input_tagalign> <output_pr1> <output_pr2> <replicate_name> <is_paired> <random_seed> <log>

set -euo pipefail

# Parse arguments
INPUT_TAGALIGN="$1"
OUTPUT_PR1="$2"
OUTPUT_PR2="$3"
REPLICATE_NAME="$4"
IS_PAIRED="$5"
RANDOM_SEED="$6"
LOG="$7"

# Create output directory
OUTPUT_DIR=$(dirname "$OUTPUT_PR1")
mkdir -p "$OUTPUT_DIR"

# Temporary files
PREFIX="$OUTPUT_DIR/$REPLICATE_NAME"
TMP_PR1="$PREFIX.00"
TMP_PR2="$PREFIX.01"

# Cleanup function
cleanup() {
    rm -f "$TMP_PR1" "$TMP_PR2"
}
trap cleanup EXIT

# Determine random seed
if [ "$RANDOM_SEED" = "0" ]; then
    random_seed=$(zcat -f "$INPUT_TAGALIGN" | wc -c)
else
    random_seed="$RANDOM_SEED"
fi

if [ "$IS_PAIRED" = "True" ]; then
    # Paired-end: keep pairs together
    nlines=$(zcat -f "$INPUT_TAGALIGN" | wc -l)
    nlines=$((nlines / 2))
    nlines=$((nlines / 2 + 1))

    # Each pair = 2 lines. We'll glue them with sed, shuffle, split, and then unglue.
    # Remove ERR trap to prevent SIGPIPE (exit 141) errors from causing hard failure:
    set +e

    zcat -f "$INPUT_TAGALIGN" | sed 'N;s/\n/\t/' | \
    shuf --random-source=<(openssl enc -aes-256-ctr -pass pass:$random_seed -nosalt </dev/zero 2>/dev/null) | \
    split -d -l $nlines - "$PREFIX." 2>> "$LOG"

    set -e

    # Restore paired-end format (convert each tabbed line back to 2 lines)
    awk -F $'\t' 'NF==12{printf "%s\t%s\t%s\t%s\t%s\t%s\n%s\t%s\t%s\t%s\t%s\t%s\n",$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' "$TMP_PR1" | gzip -nc > "$OUTPUT_PR1" 2>> "$LOG"
    awk -F $'\t' 'NF==12{printf "%s\t%s\t%s\t%s\t%s\t%s\n%s\t%s\t%s\t%s\t%s\t%s\n",$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' "$TMP_PR2" | gzip -nc > "$OUTPUT_PR2" 2>> "$LOG"
else
    # Single-end
    nlines=$(zcat -f "$INPUT_TAGALIGN" | wc -l)
    nlines=$((nlines / 2 + 1))

    # Remove ERR trap to prevent SIGPIPE (exit 141) errors from causing hard failure:
    set +e

    zcat -f "$INPUT_TAGALIGN" | \
    shuf --random-source=<(openssl enc -aes-256-ctr -pass pass:$random_seed -nosalt </dev/zero 2>/dev/null) | \
    split -d -l $nlines - "$PREFIX." 2>> "$LOG"

    set -e

    gzip -nc "$TMP_PR1" > "$OUTPUT_PR1" 2>> "$LOG"
    gzip -nc "$TMP_PR2" > "$OUTPUT_PR2" 2>> "$LOG"
fi

# Cleanup handled by trap EXIT




