rule bwa_mem_se:
	input:
		"hg19/hg19_bwa_index/genome",
		"fastq/{sample}.fastq"
	output:
		"bwa_mapped_reads/{sample}.bam"
	params: RG = "@RG\\tID:{sample}_SE\\tSM:{sample}\\tPL:Illumina"
	threads: 8
	shell:
		"""
		bwa mem -t {threads} -R '{params.RG}' {input[0]} {input[1]} | 
		samtools view -bS -@ {threads} | 
		samtools sort -o {output} -@ 8
		"""
