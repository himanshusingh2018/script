configfile: "profile/config.yaml"

rule picard_markduplicate:
	input:
		picard_clean=config['resources']['picard_clean']
	output:
		picard_mkdup=config['resources']['picard_mkdup_bam'],
		picard_metrics=config['resources']['picard_mkdup_txt']

	shell:
		"""
		picard -Xmx64g MarkDuplicates \
		--INPUT {input} \
		--OUTPUT {output[0]} \
		--METRICS_FILE {output[1]} \
		--VALIDATION_STRINGENCY LENIENT \
		--CREATE_INDEX true
		"""
