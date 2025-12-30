#!/usr/bin/env python3
"""
Prealignment script for filtering reads against contaminant genomes.
Moved from inline run: block to improve maintainability.
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
original_fq1 = str(input.fq1)
original_fq2 = str(input.fq2) if paired and input.fq2 else None
current_fq1 = original_fq1
current_fq2 = original_fq2
aligner = params.aligner.lower()

stats_path = str(output.stats)
stats_lines = []

output_fq1 = str(output.fq1)
output_fq2 = str(output.fq2)

# Handle empty prealignment list
if not prealign_list:
    shutil.copyfile(current_fq1, output_fq1)
    if paired and current_fq2:
        shutil.copyfile(current_fq2, output_fq2)
    else:
        Path(output_fq2).parent.mkdir(parents=True, exist_ok=True)
        Path(output_fq2).touch(exist_ok=True)
    total_reads = count_reads(current_fq1)
    with open(stats_path, "w") as stats_fh:
        stats_fh.write("prealignment\treads_before\treads_after\treads_filtered\tpercent_filtered\n")
        stats_fh.write(f"total\t{total_reads}\t{total_reads}\t0\t0.000000\n")
    sys.exit(0)

# Process each prealignment step
keep_set = {original_fq1, output_fq1}
if paired and original_fq2:
    keep_set.add(original_fq2)
if paired and output_fq2:
    keep_set.add(output_fq2)

for entry in prealign_list:
    name, index = parse_entry(entry)
    stage_prefix = os.path.join(prealign_dir, f"{wildcards.sample_run}.{name}.unmapped")
    prev_fq1 = current_fq1
    prev_fq2 = current_fq2
    reads_before = count_reads(prev_fq1)

    if aligner == "bowtie2":
        if paired:
            cmd = (
                "set -euo pipefail; "
                f"bowtie2 -x {quote(index)} -1 {quote(current_fq1)} -2 {quote(current_fq2)} "
                f"-p {threads} --un-conc-gz {quote(stage_prefix)} -S /dev/null "
                f"2>> {quote(log_file)}"
            )
            shell(cmd)
            current_fq1 = f"{stage_prefix}.1.gz"
            current_fq2 = f"{stage_prefix}.2.gz"
        else:
            single_target = f"{stage_prefix}.gz"
            cmd = (
                "set -euo pipefail; "
                f"bowtie2 -x {quote(index)} -U {quote(current_fq1)} "
                f"-p {threads} --un-gz {quote(single_target)} -S /dev/null "
                f"2>> {quote(log_file)}"
            )
            shell(cmd)
            current_fq1 = single_target
            current_fq2 = None
    else:  # bwa-mem2
        if paired:
            fq1_tmp = f"{stage_prefix}_1.fq"
            fq2_tmp = f"{stage_prefix}_2.fq"
            cmd = (
                "set -euo pipefail; "
                f"bwa-mem2 mem -t {threads} {quote(index)} {quote(current_fq1)} {quote(current_fq2)} "
                f"2>> {quote(log_file)} | "
                f"samtools fastq -f 4 -1 {quote(fq1_tmp)} -2 {quote(fq2_tmp)} "
                f"-0 /dev/null -s /dev/null -n - 2>> {quote(log_file)}"
            )
            shell(cmd)
            shell(f"gzip -f {quote(fq1_tmp)}")
            shell(f"gzip -f {quote(fq2_tmp)}")
            current_fq1 = f"{fq1_tmp}.gz"
            current_fq2 = f"{fq2_tmp}.gz"
        else:
            fq_tmp = f"{stage_prefix}.fq"
            cmd = (
                "set -euo pipefail; "
                f"bwa-mem2 mem -t {threads} {quote(index)} {quote(current_fq1)} "
                f"2>> {quote(log_file)} | "
                f"samtools fastq -f 4 -0 {quote(fq_tmp)} -s /dev/null -n - "
                f"2>> {quote(log_file)}"
            )
            shell(cmd)
            shell(f"gzip -f {quote(fq_tmp)}")
            current_fq1 = f"{fq_tmp}.gz"
            current_fq2 = None

    reads_after = count_reads(current_fq1)
    filtered_reads = max(reads_before - reads_after, 0)
    percent_filtered = filtered_reads / reads_before if reads_before else 0.0
    stats_lines.append((name, reads_before, reads_after, filtered_reads, percent_filtered))

    if prev_fq1 != original_fq1:
        cleanup(prev_fq1, keep_set, prealign_dir)
    if paired and prev_fq2 and prev_fq2 != original_fq2:
        cleanup(prev_fq2, keep_set, prealign_dir)

# Move final files to output
os.makedirs(os.path.dirname(output_fq1), exist_ok=True)
os.makedirs(os.path.dirname(output_fq2), exist_ok=True)

if current_fq1 != output_fq1:
    if os.path.exists(output_fq1):
        os.remove(output_fq1)
    shutil.move(current_fq1, output_fq1)

if paired and current_fq2:
    if current_fq2 != output_fq2:
        if os.path.exists(output_fq2):
            os.remove(output_fq2)
        shutil.move(current_fq2, output_fq2)
else:
    Path(output_fq2).touch(exist_ok=True)

# Write statistics
if stats_lines:
    total_before = stats_lines[0][1]
    total_after = stats_lines[-1][2]
else:
    total_before = count_reads(original_fq1)
    total_after = count_reads(output_fq1)
total_filtered = max(total_before - total_after, 0)
total_percent = total_filtered / total_before if total_before else 0.0

with open(stats_path, "w") as stats_fh:
    stats_fh.write("prealignment\treads_before\treads_after\treads_filtered\tpercent_filtered\n")
    for row in stats_lines:
        stats_fh.write(f"{row[0]}\t{row[1]}\t{row[2]}\t{row[3]}\t{row[4]:.6f}\n")
    stats_fh.write(f"total\t{total_before}\t{total_after}\t{total_filtered}\t{total_percent:.6f}\n")


