# Get runs per sample for merging
def get_runs_for_sample(sample_name):
    """Get all runs for a given sample name"""
    sample_runs = annot[annot['sample_name'] == sample_name].index.tolist()
    return sample_runs if len(sample_runs) > 0 else []

def get_units_fastqs(wildcards):
    """Get fastq files for a sample_run"""
    u = annot.loc[wildcards.sample_run]
    return [u["fq1"], u["fq2"]]

def get_quantifications(wildcards):

    if wildcards.kind=="support":
        paths = expand(os.path.join(result_path, "results", "{sample}", "peaks", "{sample}_quantification_support_counts.csv"), sample=samples_quantify)

    if wildcards.kind=="consensus":
        paths = expand(os.path.join(result_path, "results", "{sample}", "mapped", "{sample}_quantification_consensus_counts.csv"), sample=samples_quantify)

    if wildcards.kind=="promoter":
        paths = expand(os.path.join(result_path, "results", "{sample}", "mapped", "{sample}_quantification_promoter_counts.csv"), sample=samples_quantify)

    return paths