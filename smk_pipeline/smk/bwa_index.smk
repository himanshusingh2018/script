configfile: "profile/config.yaml"

rule bwa_index:
    input:
        genome=config['resources']['genome']
    output:
        bwa_index=config['resources']['bwa_index']
    shell:
        """
        bwa index -p {output.bwa_index} {input.genome} | 
        touch {output.bwa_index}
        """


