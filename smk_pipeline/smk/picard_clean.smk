configfile: "profile/config.yaml"

rule picard_clean:
	input:
		bam=config['resources']['bwa_map']
	output:
		picard_clean=config['resources']['picard_clean']
	shell:
		"""
		picard -Xmx64g \
		CleanSam --INPUT {input.bam} \
		--OUTPUT {output.picard_clean}
		"""
