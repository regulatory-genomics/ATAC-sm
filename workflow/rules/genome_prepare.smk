
rule get_genome_fai:
    input:
        genome_fasta = config["genome_fasta"]
    output:
        fai = os.path.join(result_path, "genome", "genome.fai")
    log:
        "logs/rules/get_genome_fai.log"
    conda:
        "../envs/bowtie2.yaml"
    shell:
        """
        mkdir -p $(dirname {output.fai})
        samtools faidx {input.genome_fasta}
        cp {input.genome_fasta}.fai {output.fai} 2> {log}
        """

rule get_autosomes:
    input:
        fai = os.path.join(result_path, "genome", "genome.fai")
    output:
        txt = os.path.join(result_path, "genome", "autosomes.txt"),
    log:
        "logs/get_autosomes/get_autosomes.log"
    conda:
        "../envs/bowtie2.yaml"
    shell:
        """
        # Get autosomes from FAI file
        python {workflow.basedir}/scripts/get_autosome.py {input.fai} {output.txt} 2> {log}
        """



rule gtf2tss:
    input:
        gtf = config.get("gtf", config.get("gencode_gtf", ""))
    output:
        tss = os.path.join(result_path, "genome", "tss.bed"),
    log:
        "logs/rules/gtf2tss/gtf2tss.log"
    conda:
        "../envs/macs2_homer.yaml" # Because this yaml contains perl
    shell:
        """
        mkdir -p $(dirname {output.tss})
        # Convert GTF to BED format
        perl {workflow.basedir}/scripts/gtf2tss.pl {input.gtf} > {output.tss} 2> {log}
        """