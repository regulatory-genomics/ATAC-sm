# ============================================================================
# Core Logic Simplifications
# ============================================================================

def get_runs_for_sample(sample_name):
    """Get all runs for a given sample name"""
    sample_runs = annot[annot['sample_name'] == sample_name].index.tolist()
    return sample_runs if sample_runs else []

def get_samples_passing_qc():
    """Get sample names that have at least one run passing QC (passqc=1).
    If passqc column doesn't exist, returns all samples (default: all pass QC)."""
    if 'passqc' not in annot.columns:
        # Default: all samples pass QC if column doesn't exist
        return list(samples.keys())
    
    # Get samples where at least one run has passqc == 1
    # Handle both integer and string representations
    passing_samples = set()
    for sample_name in samples.keys():
        sample_runs = get_runs_for_sample(sample_name)
        for sample_run in sample_runs:
            passqc_value = annot.loc[sample_run, 'passqc']
            # Handle int, string, or boolean values
            if passqc_value == 1 or str(passqc_value).strip() == '1':
                passing_samples.add(sample_name)
                break  # At least one run passes, include the sample
    
    return list(passing_samples) if passing_samples else []

# SIMPLIFIED: Moved nested function out for reusability
def expand_pep_path(path):
    """Expand PEP derive modifier paths like raw_data|filename"""
    if pd.isna(path) or not isinstance(path, str):
        return path
    if "|" in path:
        prefix, filename = path.split("|", 1)
        # Try config paths
        if prefix == "raw_data" and "paths" in config and "data_dir" in config["paths"]:
            return os.path.join(config["paths"]["data_dir"], filename)
        # Try PEP project config
        try:
            sources = pep_project.config.get('sample_modifiers', {}).get('derive', {}).get('sources', {})
            if prefix in sources:
                return os.path.join(sources[prefix], filename)
        except:
            pass
    # Ensure absolute path
    if path and not os.path.isabs(path):
        return os.path.abspath(path)
    return path

def get_units_fastqs(wildcards):
    """Get fastq files for a sample_run"""
    u = annot.loc[wildcards.sample_run]
    fq1 = expand_pep_path(u["R1"])
    fq2 = expand_pep_path(u["R2"])
    return [fq1, fq2 if not pd.isna(u["R2"]) else None]

# SIMPLIFIED: Consolidated file gathering logic
def _get_fastqs_for_sample(sample_name, path_func):
    """Generic helper to gather R1/R2 lists using a specific path function."""
    runs = get_runs_for_sample(sample_name)
    r1s, r2s = [], []
    for run in runs:
        p1, p2 = path_func(run)
        if p1:
            r1s.append(p1)
        if p2:
            r2s.append(p2)
    return r1s, r2s

def get_all_fastqs_for_sample(sample_name):
    """Get all R1 and R2 fastq files for a sample (raw input)."""
    # Using a lambda to adapt raw data structure to the _get_fastqs_for_sample interface
    def _raw_path_func(run):
        u = annot.loc[run]
        return str(u["R1"]) if not pd.isna(u["R1"]) else None, \
               str(u["R2"]) if not pd.isna(u["R2"]) else None
    return _get_fastqs_for_sample(sample_name, _raw_path_func)

def get_quantifications(wildcards):
    """Refactored to avoid if/else block"""
    supported_kinds = ["support", "consensus", "promoter", "TSS"]
    if wildcards.kind not in supported_kinds:
        raise ValueError(f"Unsupported kind '{wildcards.kind}'. Expected: {supported_kinds}")
    
    return expand(
        os.path.join(result_path, "downstream_res", "quantification", f"{{sample}}_quantification_{wildcards.kind}_counts.csv"),
        sample=samples.keys()
    )



# ============================================================================
# Replicate & Group Simplifications
# ============================================================================

# SIMPLIFIED: Renamed from get_reproducibility_sample and made main accessor
def get_samples_for_replicate(replicate_group):
    """Get list of sample names belonging to a replicate group"""
    if 'replicate_sample_name' not in annot.columns:
        return []
    samples_list = annot[annot['replicate_sample_name'] == replicate_group]['sample_name'].unique().tolist()
    return samples_list if len(samples_list) > 1 else []

def get_replicate_names():
    """Get list of valid replicate group names (those with >1 sample)"""
    if 'replicate_sample_name' not in annot.columns:
        return []
    all_groups = annot['replicate_sample_name'].dropna().astype(str).unique()
    return [g for g in all_groups if len(get_samples_for_replicate(g)) > 1]

# SIMPLIFIED: Generic function to avoid repeating list comprehensions
def _get_files_for_group(group_id, pattern):
    """Generic helper to get files for all samples in a group using a pattern string."""
    samples_list = get_samples_for_replicate(group_id)
    return [pattern.format(sample=s, result_path=result_path) for s in samples_list]

def get_replicate_bams(replicate_name):
    """Get BAM file paths for all samples in a replicate group"""
    return _get_files_for_group(replicate_name, "{result_path}/important_processed/bam/{sample}.filtered.bam")

def get_replicate_tagaligns(replicate_name):
    """Get tagAlign file paths for all samples in a replicate group"""
    return _get_files_for_group(replicate_name, "{result_path}/middle_files/bed/{sample}.tagAlign.gz")

def get_group_pr1_tagaligns(group_id):
    """Get all pr1 tagAlign files for replicates in a group"""
    return _get_files_for_group(group_id, "{result_path}/middle_files/replicates/{sample}/{sample}.pr1.tagAlign.gz")

def get_group_pr2_tagaligns(group_id):
    """Get all pr2 tagAlign files for replicates in a group"""
    return _get_files_for_group(group_id, "{result_path}/middle_files/replicates/{sample}/{sample}.pr2.tagAlign.gz")

# Helper function to get all combined pseudo-replicate peak files
# Returns list of peak files: [rep1-pr1_vs_rep1-pr2.narrowPeak.gz, rep2-pr1_vs_rep2-pr2.narrowPeak.gz, ...]
def get_all_replicate_peaks_pr():
    """Get list of all combined pseudo-replicate peak files for downstream analysis"""
    replicate_names = get_replicate_names()
    # Get all samples in all replicate groups
    all_samples = []
    for rep in replicate_names:
        all_samples.extend(get_samples_for_replicate(rep))
    return [
        os.path.join(result_path, "middle_files", "replicates", sample, f"{sample}-pr1_vs_{sample}-pr2.narrowPeak.gz")
        for sample in all_samples
    ]

# Get all replicate group IDs (unique groups, not individual replicates)
def get_replicate_group_ids():
    """Get unique replicate group identifiers"""
    if 'replicate_sample_name' not in annot.columns:
        return []
    # Get unique group names where each group has >1 sample
    group_names = []
    for group in annot['replicate_sample_name'].dropna().astype(str).unique():
        samples_in_group = get_samples_for_replicate(group)
        if len(samples_in_group) > 1:
            group_names.append(group)
    return group_names

# Get all pairs of replicates within a group for pairwise comparisons
def get_replicate_pairs(group_id):
    """Get all pairs (i,j) where i<j within a replicate group"""
    samples_in_group = get_samples_for_replicate(group_id)
    pairs = []
    for i in range(len(samples_in_group)):
        for j in range(i+1, len(samples_in_group)):
            pairs.append((samples_in_group[i], samples_in_group[j]))
    return pairs

# Determine if paired-end based on first sample in replicate group
def is_paired_end(replicate_name):
    """Check if replicate is paired-end (based on first sample)"""
    samples_list = get_samples_for_replicate(replicate_name)
    if not samples_list:
        return False
    # Use first sample's read_type
    first_sample = samples_list[0]
    return samples[first_sample].get("read_type", "single") == "paired"

# Validation function to check if replicate name is valid
def validate_replicate_name(wildcards):
    """Validate that the replicate wildcard corresponds to a valid replicate group"""
    valid_replicates = get_replicate_names()
    if wildcards.replicate not in valid_replicates:
        # Get all sample names for better error message
        all_sample_names = list(samples.keys())
        
        # Check if user mistakenly used a sample_name instead of replicate_sample_name
        if wildcards.replicate in all_sample_names:
            # Find the correct replicate group for this sample
            if 'replicate_sample_name' in annot.columns:
                sample_row = annot[annot['sample_name'] == wildcards.replicate]
                if not sample_row.empty:
                    correct_replicate = sample_row['replicate_sample_name'].iloc[0]
                    raise ValueError(
                        f"ERROR: '{wildcards.replicate}' is a SAMPLE name, not a replicate GROUP name!\n"
                        f"  Sample '{wildcards.replicate}' belongs to replicate group: '{correct_replicate}'\n"
                        f"  Use replicate='{correct_replicate}' instead.\n\n"
                        f"Valid replicate groups (from 'replicate_sample_name' column): {valid_replicates}\n"
                        f"  Each group contains multiple samples for reproducibility analysis."
                    )
        
        raise ValueError(
            f"ERROR: Invalid replicate name '{wildcards.replicate}'!\n"
            f"Valid replicate groups: {valid_replicates}\n"
            f"  These are values from the 'replicate_sample_name' column with >1 sample.\n"
            f"  If you see an empty list, check that:\n"
            f"    1. Your annotation has a 'replicate_sample_name' column\n"
            f"    2. At least one replicate group has >1 sample (needed for reproducibility analysis)"
        )
    return []



# ============================================================================
# Alignment Parameter Helpers
# ============================================================================

def get_bowtie2_input_string(wildcards, input):
    """Generate bowtie2 input string from FASTQ lists."""
    if samples[wildcards.sample]["read_type"] == "paired" and input.fasta_rev and len(input.fasta_rev) > 0:
        return f"-1 {','.join(input.fasta_fwd)} -2 {','.join(input.fasta_rev)}"
    return f"-U {','.join(input.fasta_fwd)}"

def get_bwa_input_string(wildcards, input):
    """Generate bash process substitution string for BWA-MEM2."""
    cmd = f"<(cat {' '.join(input.fasta_fwd)})"
    if samples[wildcards.sample]["read_type"] == "paired" and input.fasta_rev and len(input.fasta_rev) > 0:
        cmd += f" <(cat {' '.join(input.fasta_rev)})"
    return cmd

def get_filtering_flags(wildcards):
    """Generate samtools filtering flags based on read type."""
    base_flags = f"-q 30 -F {config['filtering']['sam_flag']} -L {config['refs']['whitelist']}"
    if samples[wildcards.sample]["read_type"] == "paired":
        return f"{base_flags} -f 2"
    return base_flags

def get_add_mate_tags(wildcards):
    """Get samblaster addMateTags flag if paired-end."""
    return "--addMateTags" if samples[wildcards.sample]["read_type"] == "paired" else ""

def sanitize_bwa_min_score_flag():
    """Sanitize BWA min_score parameter."""
    value = config["alignment"]["bwa"].get("min_score")
    if value is None:
        return ""
    if isinstance(value, str):
        stripped = value.strip()
        if stripped.lower() in {"", "null", "none"}:
            return ""
        return f"-T {stripped}"
    return f"-T {value}"

def get_bwa_index_path():
    """Get BWA-MEM2 index path, falling back to genome FASTA if not set."""
    index = config["alignment"]["bwa"].get("index")
    if index and index not in ("", "null", None):
        return index
    return config["refs"]["fasta"]

def get_bwa_index_input(wildcards=None):
    """Get BWA-MEM2 index file dependencies.
    Accepts wildcards parameter for Snakemake input function compatibility."""
    if not config["alignment"]["bwa"].get("index") or config["alignment"]["bwa"].get("index") in ("", "null", None):
        return multiext(config["refs"]["fasta"], ".amb", ".ann", ".bwt", ".pac", ".sa", ".0123", ".alt")
    return []

def _bwa_mem_mb(wildcards, attempt):
    """Return memory allocation for bwa alignment based on retry attempt."""
    if attempt == 1:
        return 20000  # 20G for first attempt
    elif attempt == 2:
        return 40000  # 40G for second attempt
    else:
        return 100000  # 100G for third+ attempt


# ============================================================================
# Path Management Helpers
# ============================================================================

def get_output_dir(subdir):
    """Centralized path management using Path objects."""
    return Path(result_path) / subdir

def get_trimmed_fastq_paths(sample_run):
    """Get trimmed FASTQ paths for a sample_run."""
    base_dir = get_output_dir("middle_files/trimmed")
    if annot.loc[sample_run, "read_type"] == "paired":
        return (
            str(base_dir / f"{sample_run}_1.fq.gz"),
            str(base_dir / f"{sample_run}_2.fq.gz"),
        )
    return (
        str(base_dir / f"{sample_run}.fq.gz"),
        None,
    )

def get_prealigned_fastq_paths(sample_run):
    """Get prealigned FASTQ paths for a sample_run."""
    base_dir = get_output_dir("middle_files/prealigned")
    return (
        str(base_dir / f"{sample_run}_prealigned_1.fq.gz"),
        str(base_dir / f"{sample_run}_prealigned_2.fq.gz"),
    )

def get_all_trimmed_fastqs_for_sample(sample_name):
    """Get all trimmed FASTQ files for a sample (all runs)."""
    return _get_fastqs_for_sample(sample_name, get_trimmed_fastq_paths)

def get_all_prealigned_fastqs_for_sample(sample_name):
    """Get all prealigned FASTQ files for a sample (all runs)."""
    return _get_fastqs_for_sample(sample_name, get_prealigned_fastq_paths)

def get_reads(wildcards, direction=0):
    """
    Determines input FASTQs (Trimmed vs Prealigned) based on configuration.
    direction: 0 for R1, 1 for R2.
    """
    sample = wildcards.sample
    path_func = get_prealigned_fastq_paths if has_prealignments else get_trimmed_fastq_paths
    files = _get_fastqs_for_sample(sample, path_func)
    return files[direction]