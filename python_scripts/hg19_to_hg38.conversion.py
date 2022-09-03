import os
import pandas as pd

list = [f for f in os.listdir('.') if f.endswith('.bed')]

for item in list:
    file = pd.read_csv(item,sep="\t",header=None,usecols=[0,1,2,3],names=['chr','start','end','score'])
    file.to_csv(item.split('_K562_against_Input_peaks.bed')[0]+'.hg19.bed',sep="\t",header=None,index=False)
    cmd = '/home/singh/Desktop/Software/liftOver/liftOver '+item.split('_K562_against_Input_peaks.bed')[0]+'.hg19.bed'+' /home/singh/mount/sacapus_remote/Data/liftOver_refseq/hg19ToHg38.over.chain.gz '+item.split('_K562_against_Input_peaks.bed')[0]+'.hg38.bed '+item.split('_K562_against_Input_peaks.bed')[0]+'.hg38.unlifted'
    os.system(cmd)
    print(item + ' lifted from hg19 to hg38 ...')
