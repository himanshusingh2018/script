configfile: "profile/config.yaml"

rule bowtie2_index:
    input:
        genome_fa=config['resources']['genome']
    output:
        index=config['resources']['bowtie2_index']
    conda:
        "envs/bowtie2.yaml"
    shell:
        """
        bowtie2-build {input.genome_fa} {output.index} |
        touch
        """
