import pandas as pd
import numpy as np
import os
import subprocess
import seaborn as sns
import datetime
import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager

start_time = datetime.datetime.now()
print("\nThe job started at: ",start_time,"\n")

print(
        "Following inputs are required from user:\n"
        "1: List of Induced Genes\n"
        "   ['gene']\n"
        "2: RefSeq GTF file: 'hg38.refGene.txt'\n"
     )

class Extract_induced_genes_promoter:
    def __init__(self):
        '''Read File'''
    def induced_genes(self,refseq_gtf,genelist):
        print("##### Extraction induced genes coordinates from RefSeq Database file: 'RefGene.txt' ...")
        self.genelist = genelist
        self.refseq_gtf = refseq_gtf

        def pStart(s):
            if(s['strand'] == '+'):
                return s['txStart'] - 200
            else:
                return s['txEnd'] - 50
        def pEnd(s):
            if(s['strand'] == '+'):
                return s['txStart'] + 50
            else:
                return s['txEnd'] + 200

        induced_gene_data = pd.read_csv(self.genelist,sep="\t",header=0,usecols=['gene'])

        refgene_col = ["bin","name","chrom","strand","txStart","txEnd","cdsStart","cdsEnd","exonCount","exonStarts","exonEnds","score","name2","cdsStartStat","cdsEndStat","exonFrames"]
        refgene_data = pd.read_csv(self.refseq_gtf, sep="\t", header=None,names=refgene_col)
        refgene_data = refgene_data[~refgene_data.chrom.str.contains('_')]
        refgene_data['txn_length'] = abs(refgene_data['txStart']-refgene_data['txEnd'])
        refgene_data_txn_length = refgene_data.groupby('name2', as_index=False).apply(pd.DataFrame.sort_values, 'txn_length', ascending = False).reset_index(drop=True)
        induced_refgene_data = refgene_data_txn_length[refgene_data_txn_length['name2'].isin(induced_gene_data['gene'])].copy()
        induced_refgene_data['tss'] = induced_refgene_data.groupby(['name2']).cumcount()+1
        induced_gene_tss = pd.DataFrame(induced_refgene_data, columns=['name','chrom','strand','txStart','txEnd','name2','tss'])
        induced_gene_tss['promoterStart'] = induced_gene_tss.apply(pStart, axis=1)
        induced_gene_tss['promoterEnd'] = induced_gene_tss.apply(pEnd, axis=1)
        #induced_gene_tss['promoterStart'] = np.where(induced_gene_tss.chrom == '+',induced_gene_tss['txStart']-250,induced_gene_tss['txEnd']-50)
        #induced_gene_tss['promoterEnd'] = [induced_gene_tss['txStart']-250 if induced_gene_tss.chrom == '+' else induced_gene_tss['txStart']+50 for x in induced_gene_tss]
        print(induced_gene_tss.head(4))
        induced_gene_tss.to_csv('Epromoter_IFN_k562_ENCODE_RNA_Seq_DESeq2_extracted_promoter.hg19.tsv',sep="\t",header=True,index=False)
        #print("  Induced genes tss coordinates are saved in variable: 'induced_gene_tss'\n")

x = Extract_induced_genes_promoter()
x.induced_genes('hg19.refGene.txt','RNAseq_K562_DESeq2_EdgeR_pval_and_norm_count_log2.txt')

print("The job is completed at: ",datetime.datetime.now())
print("Total time cost: ",datetime.datetime.now()-start_time)
