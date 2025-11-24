# Get runs per sample for merging
def get_runs_for_sample(sample_name):
    """Get all runs for a given sample name"""
    sample_runs = annot[annot['sample_name'] == sample_name].index.tolist()
    return sample_runs if len(sample_runs) > 0 else []

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
    return paths