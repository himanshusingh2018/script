configfile: "profile/config.yaml"

rule picard_clean:
	input:
		bam=config['resources']['bwa_map']
	output:
		picard_clean=config['resources']['picard_clean']
	
	conda: config['resources']['cwd']+'envs/picard.yaml'
	log: config['resources']['cwd']+'log/picard_clean/{sample}_log.txt'
	message: "PICARD CLEANSAM: {input.bam}"

	shell:
		"""
		picard -Xmx64g \
		CleanSam --INPUT {input.bam} \
		--OUTPUT {output.picard_clean}
		"""
