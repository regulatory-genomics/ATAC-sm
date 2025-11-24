# annotate consensus regions

# prepare configs for uropa
rule uropa_prepare:
    input:
        consensus_regions = os.path.join(result_path,"downstream_res","annotation","consensus_regions.bed"),
        gencode_template = workflow.source_path(config["annotation"]["templates"]["gencode"]),
        regulatory_template = workflow.source_path(config["annotation"]["templates"]["regulatory"]),
        gencode_gtf = config["refs"]["gencode_gtf"],
        regulatory_build_gtf = config["refs"]["regulatory_gtf"],
    output:
        gencode_config = os.path.join(result_path,"middle_files","annotation","consensus_regions_gencode.json"),
        reg_config = os.path.join(result_path,"middle_files","annotation","consensus_regions_reg.json"),
    params:
        tss_size = config["annotation"]["tss_size"],
        proximal_size_up = config["annotation"]["promoter"]["up"],
        proximal_size_dn = config["annotation"]["promoter"]["down"],
        distal_size = config["annotation"]["distal_size"],
    resources:
        mem_mb=config["resources"].get("mem_mb", 16000),
    threads: config["resources"].get("threads", 2)
    log:
        "logs/rules/annotation/uropa_prepare.log"
    run:
        ### generate gencode config
        with open(input.gencode_template) as f:
            gencode_template=Template(f.read())

        gencode_config=gencode_template.substitute({
                'TSS_flanking':'"{}"'.format(params.tss_size),
                'TSS_proximal_upstream':'"{}"'.format(params.proximal_size_up),
                'TSS_proximal_downstream':'"{}"'.format(params.proximal_size_dn),
                'distal_distance':'"{}"'.format(params.distal_size),
                'gtf_file':'"{}"'.format(input.gencode_gtf),
                'bed_file':'"{}"'.format(input.consensus_regions)
            })

        with open(output.gencode_config,'w') as out:
            out.write(gencode_config)

        ### generate reg config file
        with open(input.regulatory_template) as f:
            reg_template=Template(f.read())  

        reg_config=reg_template.substitute({
            'gtf_file':'"{}"'.format(input.regulatory_build_gtf),
            'bed_file':'"{}"'.format(input.consensus_regions)
        })

        with open(output.reg_config,'w') as out:
            out.write(reg_config)

# run uropa on consensus regions for gencode
rule uropa_gencode:
    input:
        consensus_regions = os.path.join(result_path,"downstream_res","annotation","consensus_regions.bed"),
        gencode_config = os.path.join(result_path,"middle_files","annotation","consensus_regions_gencode.json"),
    output:
        gencode_results = os.path.join(result_path,"middle_files","annotation","gencode_finalhits.txt"),
    params:
        results_dir = os.path.join(result_path,"middle_files","annotation"),
    resources:
        mem_mb=config["resources"].get("mem_mb", 16000),
    threads: 4*config["resources"].get("threads", 2)
    conda:
        "../envs/uropa.yaml",
    log:
        "logs/rules/annotation/uropa_run_gencode.log"
    shell:
        """
        uropa -p {params.results_dir}/gencode -i {input.gencode_config} -t {threads} -l {params.results_dir}/uropa.gencode.log
        """

# run uropa on consensus regions for regulatory build
rule uropa_reg:
    input:
        consensus_regions = os.path.join(result_path,"downstream_res","annotation","consensus_regions.bed"),
        reg_config = os.path.join(result_path,"middle_files","annotation","consensus_regions_reg.json"),
    output:
        reg_results = os.path.join(result_path,"middle_files","annotation","reg_finalhits.txt"),
    params:
        results_dir = os.path.join(result_path,"middle_files","annotation"),
    resources:
        mem_mb=config["resources"].get("mem_mb", 16000),
    threads: 4*config["resources"].get("threads", 2)
    conda:
        "../envs/uropa.yaml",
    log:
        "logs/rules/uropa_run_reg.log"
    shell:
        """
        uropa -p {params.results_dir}/reg -i {input.reg_config} -t {threads} -l {params.results_dir}/uropa.reg.log
        """

# peak annotation using homer
rule homer_region_annotation:
    input:
        consensus_regions = os.path.join(result_path,"downstream_res","annotation","consensus_regions.bed"),
        homer_script = os.path.join(HOMER_path,"configureHomer.pl"),
    output:
        homer_annotations = os.path.join(result_path,"middle_files","annotation","homer_annotations.tsv"),
        homer_annotations_log = os.path.join(result_path,"middle_files","annotation","homer_annotations.tsv.log"),
    params:
        homer_bin = os.path.join(HOMER_path,"bin"),
        genome = config["project"]["genome"],
    resources:
        mem_mb=config["resources"].get("mem_mb", 16000),
    threads: config["resources"].get("threads", 2)
    conda:
        "../envs/macs2_homer.yaml",
    log:
        "logs/rules/annotation/homer_region_annotation.log"
    shell:
        """
        export PATH="{params.homer_bin}:$PATH";
        
        # Check if consensus_regions file exists and is not empty
        if [ ! -f {input.consensus_regions} ] || [ ! -s {input.consensus_regions} ]; then
            # Create empty annotation files
            echo "No consensus regions found, creating empty annotation file" > {output.homer_annotations_log}
            touch {output.homer_annotations}
        else
            {params.homer_bin}/annotatePeaks.pl {input.consensus_regions} {params.genome} \
                > {output.homer_annotations} \
                2> {output.homer_annotations_log} || true;
            
            # Ensure output file exists even if command fails
            if [ ! -f {output.homer_annotations} ]; then
                touch {output.homer_annotations}
            fi
        fi
        """
        
# get gc content and region length
rule bedtools_annotation:
    input:
        consensus_regions = os.path.join(result_path,"downstream_res","annotation","consensus_regions.bed"),
        genome_fasta = config["refs"]["fasta"],
    output:
        bedtools_annotation = os.path.join(result_path, "middle_files", "annotation", "bedtools_annotation.bed"),
    resources:
        mem_mb=config["resources"].get("mem_mb", 16000),
    threads: config["resources"].get("threads", 2)
    conda:
        "../envs/pybedtools.yaml",
    log:
        "logs/rules/annotation/bedtools_annotation.log"
    shell:
        """
        bedtools nuc -fi {input.genome_fasta} -bed {input.consensus_regions} > {output.bedtools_annotation}
        """
        
# aggregate uropa and homer annotation results
rule region_annotation_aggregate:
    input:
        gencode_results = os.path.join(result_path,"middle_files","annotation","gencode_finalhits.txt"),
        reg_results = os.path.join(result_path,"middle_files","annotation","reg_finalhits.txt"),
        homer_annotations = os.path.join(result_path,"middle_files","annotation","homer_annotations.tsv"),
        bedtools_annotation = os.path.join(result_path, "middle_files", "annotation", "bedtools_annotation.bed"),
    output:
        region_annotation = os.path.join(result_path,'downstream_res','annotation',"consensus_annotation.csv"),
    resources:
        mem_mb=config["resources"].get("mem_mb", 16000),
    threads: config["resources"].get("threads", 2)
    log:
        "logs/rules/annotation/region_annotation_aggregate.log"
    run:
        import os
        import shutil
        
        # Helper function to check if file is empty
        def is_empty_file(filepath):
            return not os.path.exists(filepath) or os.path.getsize(filepath) == 0
        
        # load and format uropa gencode results
        if is_empty_file(input.gencode_results):
            gencode_characterization = pd.DataFrame(columns=['peak_chr','peak_start','peak_end','length','feat_anchor','distance','relative_location','feat_type','gene_id','gene_name','characterization'])
            gencode_characterization = gencode_characterization.set_index(pd.Index([], name='peak_id'))
        else:
            gencode_characterization = pd.read_csv(input.gencode_results, sep='\t')
            if len(gencode_characterization) == 0:
                gencode_characterization = pd.DataFrame(columns=['peak_chr','peak_start','peak_end','length','feat_anchor','distance','relative_location','feat_type','gene_id','gene_name','characterization'])
                gencode_characterization = gencode_characterization.set_index(pd.Index([], name='peak_id'))
            else:
                gencode_characterization = gencode_characterization.set_index("peak_id")
                gencode_characterization.loc[gencode_characterization['feature']=='transcript','feat_type']='transcript:'+gencode_characterization.loc[gencode_characterization['feature']=='transcript','transcript_type']
                gencode_characterization.loc[gencode_characterization['feature']=='gene','feat_type']='gene:'+gencode_characterization.loc[gencode_characterization['feature']=='gene','gene_type']
                gencode_characterization['length']=gencode_characterization['peak_end']-gencode_characterization['peak_start']
                gencode_characterization=gencode_characterization[['peak_chr','peak_start','peak_end','length','feat_anchor','distance','relative_location','feat_type','gene_id','gene_name','name']]
                gencode_characterization.columns=['chr','start','end','length','feat_anchor','distance','location','feat_type','gene_id','gene_name','characterization']
                gencode_characterization.loc[gencode_characterization['characterization'].isna(),'characterization']='NONE'
        gencode_characterization=gencode_characterization.add_prefix('gencode_')

        # load and format uropa regulatory build results
        if is_empty_file(input.reg_results):
            reg_characterization = pd.DataFrame(columns=['reg_feature','reg_feature_id'])
            reg_characterization = reg_characterization.set_index(pd.Index([], name='peak_id'))
        else:
            reg_characterization=pd.read_csv(input.reg_results,sep='\t')
            if len(reg_characterization) == 0:
                reg_characterization = pd.DataFrame(columns=['reg_feature','reg_feature_id'])
                reg_characterization = reg_characterization.set_index(pd.Index([], name='peak_id'))
            else:
                reg_characterization = reg_characterization.set_index('peak_id')[['feature','ID']]
                reg_characterization.columns=['reg_feature','reg_feature_id']
                reg_characterization.loc[reg_characterization['reg_feature'].isna(),'reg_feature']='reg_NONE'
        reg_characterization=reg_characterization.add_prefix('regulatoryBuild_')
        
        # load and format homer annotation results
        if is_empty_file(input.homer_annotations):
            homer_annotation = pd.DataFrame(columns=['Annotation','Detailed Annotation','Distance to TSS','Nearest PromoterID','Entrez ID','Nearest Unigene','Nearest Refseq','Nearest Ensembl','Gene Name','Gene Alias','Gene Description','Gene Type'])
            homer_annotation = homer_annotation.set_index(pd.Index([], name='peak_id'))
        else:
            try:
                homer_annotation = pd.read_csv(input.homer_annotations,sep='\t', index_col=0)
                if len(homer_annotation) == 0:
                    homer_annotation = pd.DataFrame(columns=['Annotation','Detailed Annotation','Distance to TSS','Nearest PromoterID','Entrez ID','Nearest Unigene','Nearest Refseq','Nearest Ensembl','Gene Name','Gene Alias','Gene Description','Gene Type'])
                    homer_annotation = homer_annotation.set_index(pd.Index([], name='peak_id'))
                else:
                    # Rename index to 'peak_id' for consistency
                    homer_annotation.index.name = 'peak_id'
                    # Check if required columns exist
                    required_cols = ['Annotation','Detailed Annotation','Distance to TSS','Nearest PromoterID','Entrez ID','Nearest Unigene','Nearest Refseq','Nearest Ensembl','Gene Name','Gene Alias','Gene Description','Gene Type']
                    available_cols = [col for col in required_cols if col in homer_annotation.columns]
                    if len(available_cols) < len(required_cols):
                        # Create missing columns with NaN
                        for col in required_cols:
                            if col not in homer_annotation.columns:
                                homer_annotation[col] = None
                    homer_annotation = homer_annotation[required_cols]
            except (pd.errors.EmptyDataError, pd.errors.ParserError):
                homer_annotation = pd.DataFrame(columns=['Annotation','Detailed Annotation','Distance to TSS','Nearest PromoterID','Entrez ID','Nearest Unigene','Nearest Refseq','Nearest Ensembl','Gene Name','Gene Alias','Gene Description','Gene Type'])
                homer_annotation = homer_annotation.set_index(pd.Index([], name='peak_id'))
        homer_annotation = homer_annotation.add_prefix('homer_')
        
        # load and format bedtools annotation results
        if is_empty_file(input.bedtools_annotation):
            bedtools_annotation = pd.DataFrame()
            bedtools_annotation = bedtools_annotation.set_index(pd.Index([], name='peak_id'))
        else:
            try:
                bedtools_annotation = pd.read_csv(input.bedtools_annotation, sep='\t', index_col = 3)
                if len(bedtools_annotation) == 0:
                    bedtools_annotation = pd.DataFrame()
                    bedtools_annotation = bedtools_annotation.set_index(pd.Index([], name='peak_id'))
                else:
                    # Rename index to 'peak_id' for consistency
                    bedtools_annotation.index.name = 'peak_id'
                    bedtools_annotation = bedtools_annotation.iloc[:,3:]
                    bedtools_annotation.columns = [col.split('_', 1)[-1].replace('at', 'AT').replace('gc', 'GC').replace('oth', 'otherBases') for col in bedtools_annotation.columns]
            except (pd.errors.EmptyDataError, pd.errors.ParserError, IndexError):
                bedtools_annotation = pd.DataFrame()
                bedtools_annotation = bedtools_annotation.set_index(pd.Index([], name='peak_id'))
        bedtools_annotation = bedtools_annotation.add_prefix('bedtools_')

        # Get common index from gencode_characterization (base)
        if len(gencode_characterization) > 0:
            base_character = gencode_characterization.copy()
        elif len(reg_characterization) > 0:
            base_character = reg_characterization.copy()
        elif len(homer_annotation) > 0:
            base_character = homer_annotation.copy()
        elif len(bedtools_annotation) > 0:
            base_character = bedtools_annotation.copy()
        else:
            # All files are empty, create empty result
            base_character = pd.DataFrame()
            base_character.index.name = 'peak_id'
        
        # join results (only if DataFrames have overlapping indices or one is empty)
        if len(base_character) > 0:
            if len(homer_annotation) > 0:
                base_character = base_character.join(homer_annotation, how='outer')
            if len(reg_characterization) > 0:
                base_character = base_character.join(reg_characterization, how='outer')
            if len(bedtools_annotation) > 0:
                base_character = base_character.join(bedtools_annotation, how='outer')
        else:
            # All empty, create empty DataFrame with all expected columns
            all_cols = list(gencode_characterization.columns) + list(reg_characterization.columns) + list(homer_annotation.columns) + list(bedtools_annotation.columns)
            base_character = pd.DataFrame(columns=all_cols)
            base_character.index.name = 'peak_id'
        
        # replace whiteapaces in colnames with underscore
        base_character.columns = base_character.columns.str.replace(' ', '_')
        
        # save final results
        base_character.to_csv(output.region_annotation, index_label='peak_id')
        
        # remove tmp folder
        middle_files_dir = os.path.dirname(input.gencode_results)
        if os.path.exists(middle_files_dir):
            shutil.rmtree(middle_files_dir)