configfile: "profile/config.yaml"

rule bam_index:
	input:
		bam=config['resources']['bwa_map']
	output:
		bai=config['resources']['bwa_map']+'.bai'
	threads: 8
	conda: config['resources']['cwd']+'envs/bwa.yaml'
	log: config['resources']['cwd']+'log/samtools_index/{sample}_log.txt'
	message: "SAMTOOLS INDEX: {input.bam}"
	shell:
		"""
		samtools index {input.bam} {output.bai} -@ {threads}
		"""
