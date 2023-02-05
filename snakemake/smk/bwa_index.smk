rule bwa_index:
	input:
		"hg19/hg19_refgene/genome.fa"
	output:
		"hg19/hg19_bwa_index/genome"
	shell:
		"""
		bwa index -p {output} {input} | 
		touch {output}
		"""

