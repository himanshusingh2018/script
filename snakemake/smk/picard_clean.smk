rule picard_clean:
	input:
		"bwa_mapped_reads/{sample}.bam"
	output:
		"picard_reads/{sample}.clean.bam"
	shell:
		"""
		picard -Xmx64g CleanSam --INPUT {input} --OUTPUT {output}
		"""
