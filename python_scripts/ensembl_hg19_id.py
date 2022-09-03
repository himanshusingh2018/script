import pandas as pd
import numpy as np

def ensembl_to_hgnc(filename):
    data = pd.read_csv(filename,sep="\t",header=None,skiprows=5,usecols=[2,8])
    data = data[data[2]=='gene']
    data['ensembl_id'] = data[8].str.split('";').str[0].str.split('"').str[1]
    data['hgnc'] = data[8].str.split('";').str[2].str.split('"').str[1]
    data = data[['ensembl_id','hgnc']]
    data.to_csv('ensembl_id_hgnc.grch37.annotation',sep="\t",header=True,index=False)

def rnaseq_annotation(rnaseq,ensembl_hgnc_annotation):
    rnaseq = pd.read_csv(rnaseq,header=0,sep="\t",usecols=['gene','RealFC'])
    rnaseq['log2FC'] = np.log2(rnaseq.RealFC)
    rnaseq = rnaseq[rnaseq.log2FC > 1]
    rnaseq.columns = ['ensembl_id','RealFC','log2FC']
    ensembl_hgnc_annotation = pd.read_csv(ensembl_hgnc_annotation, sep="\t",header=0)

    annotated_rnaseq = pd.merge(rnaseq,ensembl_hgnc_annotation,on='ensembl_id')
    annotated_rnaseq = annotated_rnaseq[['hgnc','log2FC']]
    annotated_rnaseq.columns = ['gene','log2FC']
    annotated_rnaseq.to_csv('GSE96800_PMAvsUnt.hg19.hgnc.expression',sep="\t",header=True,index=False)

ensembl_to_hgnc('/home/singh/mount/sacapus/Data/ensembl_hg19_genome_seq/Homo_sapiens.GRCh37.87.gtf')
rnaseq_annotation('GSE96800_PMAvsUnt.expression','ensembl_id_hgnc.grch37.annotation')
