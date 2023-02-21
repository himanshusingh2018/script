configfile: "profile/config.yaml"

rule bam_index:
	input:
		bam=config['resources']['bwa_map']
	output:
		bai=config['resources']['bwa_map']+'.bai'
	threads: 8
	shell:
		"""
		samtools index {input.bam} {output.bai} -@ {threads}
		"""
