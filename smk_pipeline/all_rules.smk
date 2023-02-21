configfile: "profile/config.yaml"
fq_dir = config['resources']['fq_dir']

SAMPLES = [f.split('.')[0] for f in os.listdir(fq_dir) if f.endswith('.fastq')]

rule all:
    input:
        #config['resources']['bwa_index'],
        #expand(config['resources']['bwa_map'], sample=SAMPLES),
        #expand(config['resources']['bwa_map']+'.bai', sample=SAMPLES),
        #expand(config['resources']['picard_clean'], sample=SAMPLES),
        #expand(config['resources']['picard_mkdup_bam'], sample=SAMPLES),
        expand(config['resources']['picard_mkdup_txt'], sample=SAMPLES)
        #"vcf/somatic.vcf.gz"
    output:
        "all_complete.txt"
    shell:
        "touch all_complete.txt"

include: "smk/bwa_index.smk"
include: "smk/bwa_mem_se.smk"
include: "smk/samtools_index.smk"
include: "smk/picard_clean.smk"
include: "smk/picard_markduplicate.smk"
#include: "smk/mutect2_tumor_normal_somatic_mut.smk"


