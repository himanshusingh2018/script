configfile: "profile/config.yaml"

rule bwa_mem_se:
	input:
		genome=config['resources']['bwa_index'],
		fq=config['resources']['fq_dir']+config['resources']['fq1']

	output:
		bwa_map=config['resources']['bwa_map']

	params: RG = config['resources']['RG']
	threads: 8
	shell:
		"""
		bwa mem -t {threads} -R '{params.RG}' {input.genome} {input.fq} | 
		samtools view -bS -@ {threads} | 
		samtools sort -o {output.bwa_map} -@ 8
		"""
