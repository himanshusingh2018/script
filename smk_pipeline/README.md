#smk_pipeline
1. envs: individual conda environments' yaml files
	bowtie2.yaml: bowtie
	bwa.yaml: bwa and samtools
	picard.yaml: picard

2. profile:
	config.yaml: all resources path and other parameters
3. smk
	bowtie2_index.smk: bowtie2 genome index 
	bwa_index.smk: BWA genome index
	bwa_mem_pe.smk: BWA FastQ Paired End Alignment, output sorted bam
	bwa_mem_se.smk: BWA FastQ Single End Alignment, output sorted bam
	mutect2_tumor_normal_somatic_mut.smk: GATK4 Mutect2 Variant Calling, Multiple Tumor & Single Normal, output: combined vcd.gz file
	picard_clean.smk: Picard Clean bam file
	picard_markduplicate.smk: Picard Markduplicate picard clean bam file
	samtools_index.smk: index bam file

4. all_rules.smk: call all the smk/scripts.smk