#!/usr/bin/env python3
"""
Prealignment script for filtering reads against contaminant genomes.
Moved from inline run: block to improve maintainability.
Now supports multiple input FASTQ files (all runs for a sample).
"""
import os
import sys
import shutil
import gzip
from pathlib import Path
from shlex import quote
from snakemake.shell import shell  # Standard import for shell function

# Access Snakemake objects
input = snakemake.input
output = snakemake.output
params = snakemake.params
wildcards = snakemake.wildcards
threads = snakemake.threads
log = snakemake.log

def count_reads(path):
    """Count reads in a FASTQ file."""
    if not path or not os.path.exists(path):
        return 0
    opener = gzip.open if path.endswith(".gz") else open
    total_lines = 0
    with opener(path, "rt") as handle:
        for _ in handle:
            total_lines += 1
    return total_lines // 4

def count_reads_multiple(paths):
    """Count total reads across multiple FASTQ files."""
    total = 0
    for path in paths:
        if path and os.path.exists(path):
            total += count_reads(path)
    return total

def parse_entry(entry):
    """Parse prealignment entry (dict or string)."""
    if isinstance(entry, dict):
        idx = entry.get("path")
        if idx is None:
            raise ValueError("Prealignment dictionary entries require a 'path' key.")
        name = entry.get("name") or Path(idx).stem
        return name, idx
    entry = str(entry)
    if "=" in entry:
        name, idx = entry.split("=", 1)
        return name.strip(), idx.strip()
    return Path(entry).stem, entry

def cleanup(path, keep_set, prealign_dir):
    """Remove temporary file if not in keep set."""
    if not path or not isinstance(path, str):
        return
    if path.startswith(prealign_dir) and os.path.exists(path) and path not in keep_set:
        try:
            os.remove(path)
        except:
            pass

# Main execution
prealign_list = params.prealignments
log_file = str(log)
os.makedirs(os.path.dirname(log_file), exist_ok=True)
with open(log_file, "w"):
    pass

prealign_dir = str(params.prealign_dir)
os.makedirs(prealign_dir, exist_ok=True)

paired = params.is_paired
aligner = params.aligner.lower()

stats_path = str(output.stats)
stats_lines = []

output_fq1 = str(output.fq1)
output_fq2 = str(output.fq2)

# Handle multiple input files (all runs for a sample)
input_fq1_list = input.fq1 if isinstance(input.fq1, list) else [input.fq1]
input_fq2_list = input.fq2 if isinstance(input.fq2, list) else [input.fq2]

# Filter out empty/None values
input_fq1_list = [str(f) for f in input_fq1_list if f and os.path.exists(str(f))]
input_fq2_list = [str(f) for f in input_fq2_list if f and os.path.exists(str(f))]

# Store original inputs (can be multiple files - aligners will handle them)
original_fq1_list = input_fq1_list
original_fq2_list = input_fq2_list if paired else []

# Current files start as the original inputs
current_fq1_list = original_fq1_list
current_fq2_list = original_fq2_list

# Handle empty prealignment list
if not prealign_list:
    # Concatenate multiple input files directly to output
    os.makedirs(os.path.dirname(output_fq1), exist_ok=True)
    os.makedirs(os.path.dirname(output_fq2), exist_ok=True)
    
    if len(current_fq1_list) == 1:
        shutil.copyfile(current_fq1_list[0], output_fq1)
    else:
        shell(f"cat {' '.join(quote(f) for f in current_fq1_list)} > {quote(output_fq1)}")
    
    if paired and current_fq2_list:
        if len(current_fq2_list) == 1:
            shutil.copyfile(current_fq2_list[0], output_fq2)
        else:
            shell(f"cat {' '.join(quote(f) for f in current_fq2_list)} > {quote(output_fq2)}")
    else:
        Path(output_fq2).touch(exist_ok=True)
    
    total_reads = count_reads_multiple(current_fq1_list)
    with open(stats_path, "w") as stats_fh:
        stats_fh.write("prealignment\treads_before\treads_after\treads_filtered\tpercent_filtered\n")
        stats_fh.write(f"total\t{total_reads}\t{total_reads}\t0\t0.000000\n")
    sys.exit(0)

# Process each prealignment step
keep_set = set(original_fq1_list + original_fq2_list + [output_fq1, output_fq2])

for entry in prealign_list:
    name, index = parse_entry(entry)
    stage_prefix = os.path.join(prealign_dir, f"{wildcards.sample}.{name}.unmapped")
    prev_fq1_list = current_fq1_list
    prev_fq2_list = current_fq2_list
    reads_before = count_reads_multiple(prev_fq1_list)

    if aligner == "bowtie2":
        if paired:
            # bowtie2 supports comma-separated lists for multiple files
            fq1_input = ','.join(current_fq1_list)
            fq2_input = ','.join(current_fq2_list)
            cmd = (
                "set -euo pipefail; "
                f"bowtie2 -x {quote(index)} -1 {fq1_input} -2 {fq2_input} "
                f"-p {threads} --un-conc-gz {quote(stage_prefix)} -S /dev/null "
                f"2>> {quote(log_file)}"
            )
            shell(cmd)
            current_fq1_list = [f"{stage_prefix}.1.gz"]
            current_fq2_list = [f"{stage_prefix}.2.gz"]
        else:
            # bowtie2 supports comma-separated lists for multiple files
            fq_input = ','.join(current_fq1_list)
            single_target = f"{stage_prefix}.gz"
            cmd = (
                "set -euo pipefail; "
                f"bowtie2 -x {quote(index)} -U {fq_input} "
                f"-p {threads} --un-gz {quote(single_target)} -S /dev/null "
                f"2>> {quote(log_file)}"
            )
            shell(cmd)
            current_fq1_list = [single_target]
            current_fq2_list = []
    else:  # bwa-mem2
        if paired:
            # bwa-mem2 uses process substitution with cat for multiple files
            fq1_input = f"<(cat {' '.join(quote(f) for f in current_fq1_list)})"
            fq2_input = f"<(cat {' '.join(quote(f) for f in current_fq2_list)})"
            fq1_tmp = f"{stage_prefix}_1.fq"
            fq2_tmp = f"{stage_prefix}_2.fq"
            cmd = (
                "set -euo pipefail; "
                f"bwa-mem2 mem -t {threads} {quote(index)} {fq1_input} {fq2_input} "
                f"2>> {quote(log_file)} | "
                f"samtools fastq -f 4 -1 {quote(fq1_tmp)} -2 {quote(fq2_tmp)} "
                f"-0 /dev/null -s /dev/null -n - 2>> {quote(log_file)}"
            )
            shell(cmd)
            shell(f"gzip -f {quote(fq1_tmp)}")
            shell(f"gzip -f {quote(fq2_tmp)}")
            current_fq1_list = [f"{fq1_tmp}.gz"]
            current_fq2_list = [f"{fq2_tmp}.gz"]
        else:
            # bwa-mem2 uses process substitution with cat for multiple files
            fq_input = f"<(cat {' '.join(quote(f) for f in current_fq1_list)})"
            fq_tmp = f"{stage_prefix}.fq"
            cmd = (
                "set -euo pipefail; "
                f"bwa-mem2 mem -t {threads} {quote(index)} {fq_input} "
                f"2>> {quote(log_file)} | "
                f"samtools fastq -f 4 -0 {quote(fq_tmp)} -s /dev/null -n - "
                f"2>> {quote(log_file)}"
            )
            shell(cmd)
            shell(f"gzip -f {quote(fq_tmp)}")
            current_fq1_list = [f"{fq_tmp}.gz"]
            current_fq2_list = []

    reads_after = count_reads_multiple(current_fq1_list)
    filtered_reads = max(reads_before - reads_after, 0)
    percent_filtered = filtered_reads / reads_before if reads_before else 0.0
    stats_lines.append((name, reads_before, reads_after, filtered_reads, percent_filtered))

    # Clean up previous stage files (but not original inputs)
    for f in prev_fq1_list:
        if f not in original_fq1_list:
            cleanup(f, keep_set, prealign_dir)
    for f in prev_fq2_list:
        if f not in original_fq2_list:
            cleanup(f, keep_set, prealign_dir)

# Move final files to output
os.makedirs(os.path.dirname(output_fq1), exist_ok=True)
os.makedirs(os.path.dirname(output_fq2), exist_ok=True)

# Final output should be a single file (result of last prealignment stage)
if len(current_fq1_list) == 1 and current_fq1_list[0] != output_fq1:
    if os.path.exists(output_fq1):
        os.remove(output_fq1)
    shutil.move(current_fq1_list[0], output_fq1)
elif len(current_fq1_list) > 1:
    # Shouldn't happen after prealignment stages, but handle it
    shell(f"cat {' '.join(quote(f) for f in current_fq1_list)} > {quote(output_fq1)}")
    for f in current_fq1_list:
        if f not in original_fq1_list:
            cleanup(f, keep_set, prealign_dir)

if paired and current_fq2_list:
    if len(current_fq2_list) == 1 and current_fq2_list[0] != output_fq2:
        if os.path.exists(output_fq2):
            os.remove(output_fq2)
        shutil.move(current_fq2_list[0], output_fq2)
    elif len(current_fq2_list) > 1:
        # Shouldn't happen after prealignment stages, but handle it
        shell(f"cat {' '.join(quote(f) for f in current_fq2_list)} > {quote(output_fq2)}")
        for f in current_fq2_list:
            if f not in original_fq2_list:
                cleanup(f, keep_set, prealign_dir)
else:
    Path(output_fq2).touch(exist_ok=True)

# Write statistics
if stats_lines:
    total_before = stats_lines[0][1]
    total_after = stats_lines[-1][2]
else:
    total_before = count_reads_multiple(original_fq1_list)
    total_after = count_reads(output_fq1)
total_filtered = max(total_before - total_after, 0)
total_percent = total_filtered / total_before if total_before else 0.0

with open(stats_path, "w") as stats_fh:
    stats_fh.write("prealignment\treads_before\treads_after\treads_filtered\tpercent_filtered\n")
    for row in stats_lines:
        stats_fh.write(f"{row[0]}\t{row[1]}\t{row[2]}\t{row[3]}\t{row[4]:.6f}\n")
    stats_fh.write(f"total\t{total_before}\t{total_after}\t{total_filtered}\t{total_percent:.6f}\n")

