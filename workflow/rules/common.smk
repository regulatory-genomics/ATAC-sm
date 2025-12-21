# Get runs per sample for merging
def get_runs_for_sample(sample_name):
    """Get all runs for a given sample name"""
    sample_runs = annot[annot['sample_name'] == sample_name].index.tolist()
    return sample_runs if len(sample_runs) > 0 else []

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


def get_reproducibility_sample(sample_rep):
    sample = annot['sample_name'][annot['replicate_sample_name'] == sample_rep].unique().tolist()
    return sample if len(sample) > 1 else []

def get_units_fastqs(wildcards):
    """Get fastq files for a sample_run"""
    u = annot.loc[wildcards.sample_run]
    fq1 = u["R1"]
    fq2 = u["R2"]
    
    # Handle PEP derive modifier expansion (raw_data| prefix)
    # PEP should expand this automatically, but if not, expand it manually
    def expand_pep_path(path):
        """Expand PEP derive modifier paths like raw_data|filename"""
        if pd.isna(path) or not isinstance(path, str):
            return path
        if "|" in path:
            prefix, filename = path.split("|", 1)
            # Try config paths first (already merged from PEP)
            if prefix == "raw_data" and "paths" in config and "data_dir" in config["paths"]:
                return os.path.join(config["paths"]["data_dir"], filename)
            # Fallback: try to get from PEP project config directly
            try:
                if 'sample_modifiers' in pep_project.config:
                    sources = pep_project.config.get('sample_modifiers', {}).get('derive', {}).get('sources', {})
                    if prefix in sources:
                        return os.path.join(sources[prefix], filename)
            except:
                pass
        # Ensure absolute path
        if path and not os.path.isabs(path):
            return os.path.abspath(path)
        return path
    
    fq1 = expand_pep_path(fq1)
    fq2 = expand_pep_path(fq2)
    
    if pd.isna(fq2):
        fq2 = None
    
    return [fq1, fq2]

def get_all_fastqs_for_sample(sample_name):
    """Get all R1 and R2 fastq files for a sample (all runs combined)"""
    sample_runs = get_runs_for_sample(sample_name)
    r1_files = []
    r2_files = []
    for sr in sample_runs:
        u = annot.loc[sr]
        r1 = u["R1"]
        if not pd.isna(r1):
            r1_files.append(str(r1))
        r2 = u["R2"]
        if not pd.isna(r2):
            r2_files.append(str(r2))
    return r1_files, r2_files

def get_quantifications(wildcards):
    if wildcards.kind=="support":
        paths = expand(os.path.join(result_path, "downstream_res", "quantification", "{sample}_quantification_support_counts.csv"), sample=samples.keys())
    elif wildcards.kind=="consensus":
        paths = expand(os.path.join(result_path, "downstream_res", "quantification", "{sample}_quantification_consensus_counts.csv"), sample=samples.keys())
    elif wildcards.kind=="promoter":
        paths = expand(os.path.join(result_path, "downstream_res", "quantification", "{sample}_quantification_promoter_counts.csv"), sample=samples.keys())
    elif wildcards.kind=="TSS":
        paths = expand(os.path.join(result_path, "downstream_res", "quantification", "{sample}_quantification_TSS_counts.csv"), sample=samples.keys())
    else:
        raise ValueError(f"Unsupported quantification kind '{wildcards.kind}'. Expected one of support, consensus, promoter, TSS.")
    return paths



# ============================================================================
# Helper Functions
# ============================================================================

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
    # Get all samples in all replicate groups
    all_samples = []
    for rep in replicate_names:
        all_samples.extend(get_reproducibility_sample(rep))
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
        samples_in_group = get_reproducibility_sample(group)
        if len(samples_in_group) > 1:
            group_names.append(group)
    return group_names

# Get all pairs of replicates within a group for pairwise comparisons
def get_replicate_pairs(group_id):
    """Get all pairs (i,j) where i<j within a replicate group"""
    samples_in_group = get_reproducibility_sample(group_id)
    pairs = []
    for i in range(len(samples_in_group)):
        for j in range(i+1, len(samples_in_group)):
            pairs.append((samples_in_group[i], samples_in_group[j]))
    return pairs

# Get all pr1 tagAlign files for a group
def get_group_pr1_tagaligns(group_id):
    """Get all pr1 tagAlign files for replicates in a group"""
    samples_in_group = get_reproducibility_sample(group_id)
    return [
        os.path.join(result_path, "middle_files", "replicates", sample, f"{sample}.pr1.tagAlign.gz")
        for sample in samples_in_group
    ]

# Get all pr2 tagAlign files for a group
def get_group_pr2_tagaligns(group_id):
    """Get all pr2 tagAlign files for replicates in a group"""
    samples_in_group = get_reproducibility_sample(group_id)
    return [
        os.path.join(result_path, "middle_files", "replicates", sample, f"{sample}.pr2.tagAlign.gz")
        for sample in samples_in_group
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
