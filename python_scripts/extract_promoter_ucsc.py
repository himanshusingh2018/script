#Extract promoter regions mouse and human from UCSC browser

#Genome Annotation UCSC
#mm9_gtf = "wget -c -O mm9.refGene.txt.gz http://hgdownload.soe.ucsc.edu/goldenPath/mm9/database/refGene.txt.gz"
#mm10_gtf = "wget -c -O mm10.refGene.txt.gz http://hgdownload.soe.ucsc.edu/goldenPath/mm10/database/refGene.txt.gz"
#hg19_gtf = "wget -c -O hg19.refGene.txt.gz http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz"
#hg38_gtf = "wget -c -O hg38.refGene.txt.gz http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/refGene.txt.gz"

mm9_gtf = "http://hgdownload.soe.ucsc.edu/goldenPath/mm9/database/refGene.txt.gz"
mm10_gtf = "http://hgdownload.soe.ucsc.edu/goldenPath/mm10/database/refGene.txt.gz"
hg19_gtf = "http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz"
hg38_gtf = "http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/refGene.txt.gz"


import pandas as pd
import numpy as np

def extract_promoter(genome_gtf,outPut):
    refgene_col = ["bin","name","chrom","strand","txStart","txEnd","cdsStart","cdsEnd","exonCount","exonStarts","exonEnds","score","name2","cdsStartStat","cdsEndStat","exonFrames"]
    refgene_data = pd.read_csv(genome_gtf,sep="\t",header=None,names=refgene_col,usecols=["name","chrom","strand","txStart","txEnd","name2"])
    refgene_data['startPromoter'] = np.where(refgene_data['strand'] == '+', refgene_data['txStart']-1000, refgene_data['txStart']+1000)
    refgene_data['endPromoter'] = np.where(refgene_data['strand'] == '+', refgene_data['txEnd']+1000, refgene_data['txStart']-1000)
    refgene_data.to_csv(outPut,sep="\t",header=True, index=False)
    print(outPut + ' is generated...')

extract_promoter(hg38_gtf,'hg38.all.promoters.txt')
extract_promoter(hg19_gtf,'hg19.all.promoters.txt')