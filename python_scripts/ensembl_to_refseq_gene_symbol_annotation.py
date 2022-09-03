import pandas as pd
import os

def ensembl_refseq_gene_symbol(ensembl_gtf_file_weblink):
    df = pd.read_csv(ensembl_gtf_file_weblink,skiprows=5,header=None,sep="\t")
    df = df[8].str.split(';', n = 5, expand = True)
    df.dropna(subset = [4], inplace=True)
    df = df[df[4].str.contains("gene_name")]
    df = df[[0,4]]
    df[0] = df[0].str.replace('gene_id "','')
    df[0] = df[0].str.replace('"','')
    df[4] = df[4].str.replace('gene_name "','')
    df[4] = df[4].str.replace('"','')
    df.columns = ['ENSEMBL_ID','Gene_Symbol']
    df.drop_duplicates(subset ="ENSEMBL_ID", keep = False, inplace = True)

    df.to_csv(ensembl_gtf_file_weblink.split('/')[-1].split('.')[1]+'.ensembl_to_refseq_gene_symbol.tsv',sep="\t",header=True,index=False)
    print("ENSEMBL ID to RefSeq Gene Symbol Conversion is completed...:\t\n"+os.getcwd()+'/'+ensembl_gtf_file_weblink.split('/')[-1].split('.')[1]+".ensembl_to_refseq_gene_symbol.tsv")

ensembl_refseq_gene_symbol('ftp://ftp.ensembl.org/pub/release-99/gtf/mus_musculus/Mus_musculus.GRCm38.99.gtf.gz')
