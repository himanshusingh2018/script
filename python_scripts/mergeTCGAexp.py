import pandas as pd
import glob

#delim_whitespace=True

all_exp = pd.read_csv('data/0a2fdcec-904b-4015-91bb-71e689af12a1/b033c85a-2395-4ced-a84d-ac1fbce674bf.htseq.counts.gz',sep="\t",header=None,names=['ENSG'], usecols=['ENSG'])
all_exp = all_exp[all_exp['ENSG'].str.contains('ENS') == True]


for f in glob.glob('data/*/*.gz'):
    print(f)
    tmp = pd.read_csv(f, sep="\t", header=None, names = ['ENSG',str(  f.split('/')[-1].split('.htseq')[0] )], dtype='unicode', error_bad_lines=False,engine='python')
    tmp = tmp[tmp['ENSG'].str.contains('ENSG') == True]
    all_exp = pd.merge(all_exp, tmp, on = 'ENSG', how = 'outer')
    #print( f.split('/')[-1].split('.htseq')[0] )
    #print(all_exp.head())
    
all_exp.to_csv('tcga_exp.htseq.12.23.2021.txt',sep="\t",header=True,index=False)
