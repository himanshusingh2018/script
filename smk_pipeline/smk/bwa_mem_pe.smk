configfile: "profile/config.yaml"

rule bwa_mem_pe:
	input:
		genome=config['resources']['bwa_index'],
		fq1=config['resources']['fq_dir']+config['resources']['fq1'],
		fq2=config['resources']['fq_dir']+config['resources']['fq2']

	output:
		bwa_map=config['resources']['bwa_map']

	params: RG = config['resources']['RG']
	threads: 8
	conda: config['resources']['cwd']+'envs/bwa.yaml'
	log: config['resources']['cwd']+'log/bwa_mem_se/{sample}_log.txt'
	message: "BWA MEM PE: {input.fq1} and {input.fq2}"

	shell:
		"""
		bwa mem -t {threads} -R '{params.RG}' {input.genome} {input.fq1} {input.fq2} | 
		samtools view -bS -@ {threads} | 
		samtools sort -o {output.bwa_map} -@ 8
		"""
