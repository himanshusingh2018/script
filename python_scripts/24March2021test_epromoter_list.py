import pandas as pd
import numpy as np
import sys

def read_cluster(filename):
    cluster = pd.read_csv(filename,sep="\t",header=0)
    cluster = cluster[(cluster['number of promoter binding TFs in cluster']>0)&(cluster['number of genes per cluster']>cluster['number of promoter binding TFs in cluster'])]
    cluster = cluster[['cluster number','genes in cluster']]
    #print(cluster);exit(0)
    return(cluster)
#print(read_cluster('cJun_gene_promoter_cluster_tf_binding.table'))
#exit(0)
def read_gene_tf_peak(filename):
    peak = pd.read_csv(filename,sep="\t",header=0)
    peak = peak[peak['distance'] < 1001]
    #peak = peak[peak['tf_gene_tss_distance'] < 1001]
    #print(peak.head(2));exit(0)
    return(peak)
#print(read_gene_tf_peak('cJun_gene_tss_peak.table'))
#exit(0)

def read_gene_tss_data(filename):
	gene_tss_table = pd.read_csv(filename,sep="\t",header=0)
	#print(gene_tss_table)
	gene_tss_table['txEnd'],gene_tss_table['txStart'] = np.where(gene_tss_table.strand == '-', [gene_tss_table.txStart,gene_tss_table.txEnd],[gene_tss_table.txEnd,gene_tss_table.txStart])
	gene_tss_table = gene_tss_table[['chrom','txStart','name2','tss']]
	gene_tss_table.columns = ['chr','start','gene','tss_no']
	#print(gene_tss_table.head(4))
	return(gene_tss_table)
#read_gene_tss_data('Hogan_2017_Il1b_induced_gene_tss')
#exit(0)

def epromoter_cluster(df1,df2,df3,condition,tfname,resources,outputfile):
    df = pd.DataFrame(columns=['genes in cluster', 'epromoter gene'])
    for row in df1['genes in cluster']:
        gene_list = row.split(', ')
        for item in gene_list:
            if item in df2['gene'].to_list():
                df = df.append({'genes in cluster' : row , 'epromoter gene' : item} , ignore_index=True)
    #print(df3);exit(0)
    final_df = pd.merge(df,df2,how = 'left',left_on='epromoter gene',right_on='gene')
    final_df = final_df.loc[final_df.groupby('epromoter gene')['distance'].idxmin()].reset_index(drop=True)
    #final_df = final_df.loc[final_df.groupby('epromoter gene')['tf_gene_tss_distance'].idxmin()].reset_index(drop=True)
    final_df = final_df[['genes in cluster','epromoter gene','tss_no','distance']]
    #final_df = final_df[['genes in cluster','epromoter gene','tss_no','tf_gene_tss_distance']]
    #final_df['tss_no'] = 'tss'+final_df['tss_no'].astype('str')

    final = pd.merge(df3,final_df,left_on=['gene','tss_no'],right_on=['epromoter gene','tss_no'])
    final['epromoter coordinates'] = final[['chr','start']].apply(lambda x : '{}: {}'.format(x[0],x[1]), axis=1)
    final = final[['epromoter gene','epromoter coordinates','genes in cluster','distance']]
    #final = final[['epromoter gene','epromoter coordinates','genes in cluster','tf_gene_tss_distance']]

    final = final.groupby(['genes in cluster'], as_index=False).agg(lambda x: ', '.join(set(x.dropna())))
    final['condition'] = condition
    final['TF'] = tfname
    final['Resources'] = resources
    final = final[['epromoter gene','epromoter coordinates','genes in cluster','condition','TF','Resources']]
    #print(final.head(2))
    final.to_csv(outputfile,sep="\t",header=True,index=False)
    print(outputfile, 'is generated...')

gene_promoter_cluster_tf_binding_table = str(sys.argv[1])
gene_peak_table = str(sys.argv[2])
gene_tss_data_table = str(sys.argv[3])
condition = str(sys.argv[4])
tf = str(sys.argv[5])
resource = str(sys.argv[6])
outputfile = '..\\output1\\'+str(sys.argv[7])

epromoter_cluster(read_cluster(gene_promoter_cluster_tf_binding_table),
                  read_gene_tf_peak(gene_peak_table),
                  read_gene_tss_data(gene_tss_data_table),
                  condition,
                  tf,
                  resource,
                  outputfile)

#RUN COMMAND
#python ..\\24March2021test_epromoter_list.py p65.Il1b_gene_promoter_cluster_tf_binding.table p65.Il1b_gene_tss_peak.table Hogan_2017_Il1b_induced_gene_tss Il1b p65 Hogan_2017 output.tsv