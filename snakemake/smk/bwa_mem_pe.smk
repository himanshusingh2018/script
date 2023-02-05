rule bwa_mem_pe:
	input:
		"fq/{sample}1.fastq",
		"fq/{sample}2.fastq"
	output:
		"bwa_mapped_reads/{sample}.bam"
	params: RG = "@RG\\tID:{sample}_SE\\tSM:{sample}\\tPL:Illumina",
		ref = "hg19/hg19_bwa_index/genome"
	threads: 8
	shell:
		"""
		bwa mem -t {threads} -R '{params.RG}' {params.ref} {input[0]} {input[1]} | 
		samtools view -bS -@ {threads} | 
		samtools sort -o {output} -@ {threads}
		"""
