configfile: "profile/config.yaml"

rule bwa_index:
    input:
        genome=config['resources']['genome']
    output:
        bwa_index=config['resources']['bwa_index']
    conda: config['resources']['cwd']+'envs/bwa.yaml'
    log: config['resources']['cwd']+"logs/bwa_index/logs.txt"
    message: "BWA INDEX: {input.genome}"
    shell:
        """
        bwa index -p {output.bwa_index} {input.genome} | 
        touch {output.bwa_index}
        """


