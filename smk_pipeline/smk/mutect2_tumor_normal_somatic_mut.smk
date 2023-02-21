rule mutect2_tumor_normal_somatic_mut:
	input:
		"picard_reads/N1.mkdup.bam",
		"picard_reads/T1.mkdup.bam",
		"picard_reads/T2.mkdup.bam"
	output:
		"vcf/somatic.vcf.gz"
	params:
		"hg19/hg19_refgene/genome.fa"
	shell:
		"""
		gatk Mutect2 --java-options -Xmx80g \
		-R {params[0]} \
		-I {input[0]} \
		-I {input[1]} \
		-I {input[2]} \
		-normal N1 \
		-O {output} \
		--native-pair-hmm-threads 80 
		"""
