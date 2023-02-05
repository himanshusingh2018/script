SAMPLES = [f.split('.')[0] for f in os.listdir('fastq') if f.endswith('.fastq')]

rule all:
    input:
        "hg19/hg19_bwa_index/genome",
        expand("bwa_mapped_reads/{sample}.bam", sample=SAMPLES),
        expand("bwa_mapped_reads/{sample}.bam.bai", sample=SAMPLES),
        expand("picard_reads/{sample}.clean.bam", sample=SAMPLES),
	expand("picard_reads/{sample}.mkdup.bam", sample=SAMPLES),
        expand("picard_reads/{sample}.dup.metrics.txt", sample=SAMPLES),
        "vcf/somatic.vcf.gz"
    output:
        "all_complete.txt"
    shell:
        "touch all_complete.txt"

include: "smk/bwa_index.smk"
include: "smk/bwa_mem_se.smk"
include: "smk/samtools_index.smk"
include: "smk/picard_clean.smk"
include: "smk/picard_markduplicate.smk"
include: "smk/mutect2_tumor_normal_somatic_mut.smk"

