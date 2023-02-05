rule bam_index:
	input:
		"bwa_mapped_reads/{sample}.bam"
	output:
		"bwa_mapped_reads/{sample}.bam.bai"
	threads: 8
	shell:
		"""
		samtools index {input} {output} -@ {threads}
		"""
