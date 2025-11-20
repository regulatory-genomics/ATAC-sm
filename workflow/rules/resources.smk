# install homer from here http://homer.ucsd.edu/homer/
rule install_homer:
    output:
        homer_script = os.path.join(HOMER_path,"configureHomer.pl"),
    params:
        homer_dir = HOMER_path,
        homer_url = "http://homer.ucsd.edu/homer/configureHomer.pl",
        genome = config["project"]["genome"],
    resources:
        mem_mb=config["resources"].get("mem_mb", 8000),
    threads: config["resources"].get("threads", 1)
    conda:
        "../envs/macs2_homer.yaml",
    log:
        "logs/rules/install_homer.log"
    shell:
        """
        echo "start HOMER installation (~10 minutes)"
        
        mkdir -p {params.homer_dir}
        cd {params.homer_dir}
        wget --directory-prefix={params.homer_dir} {params.homer_url}
        perl configureHomer.pl -install
        perl configureHomer.pl -install {params.genome}
        
        echo "finished HOMER installation"
        """