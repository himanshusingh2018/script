import pandas as pd
import numpy as np
import datetime
'''
Input files:
	1: file information about binding sites (chr,start,end): 'file.wig'
	2: file information with SNP information (chr, rsid,position): snp.bed
Output:
	SNPs overlapping with binding sites
	
'''

start_time = datetime.datetime.now()
print("\nThe job started at: ",datetime.datetime.now(),"\n")

def read_file(filename):#function read file
    file = pd.read_csv(filename,sep="\t",header=None)
    return file

file = read_file('file.wig')#reading file.wig
file.columns = ['chr','start','end','signal'] #assign column names file.wig
snp = read_file('snp.bed')  #read snp.bed
snp.columns = ['chr','rsid','position'] #assign columns of snp.bed

df = pd.merge(file,snp,on='chr') #merge two df on the basis of chr common column
print(df.shape)#print dimension of df
print(df.columns,"\n") #print df columns name
df.start,df.end = np.where(df.start > df.end, [df.end,df.start], [df.start,df.end]) #swap columns if start > end
df['enhancer'] = df.position.between(df.start,df.end,inclusive=True).map({True:'Yes',False:'No'}) #add new column if position is between start and end
df2=df[df.enhancer=='Yes'] #extract rows with overlapping enhancer region
df2 = df2[['chr','rsid','position','enhancer']] #extract specific columns
df2.to_csv('snp_overlapping_binding.sites',sep="\t",index=False,header=True)
print('The output file generated:\n\tsnp_overlapping_binding.sites')
print("The job is finished by: ",datetime.datetime.now())
print("Total time cost: ",datetime.datetime.now()-start_time,"\n")
