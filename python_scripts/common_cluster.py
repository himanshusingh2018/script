import pandas as pd

def read_file(filename):
    file = pd.read_csv(filename,sep="\t",header=0,usecols=['cluster number','genes in cluster'])
    file['resource'] = filename
    return(file)
    #print(file.head(4))

def common_cluster_two_dataset(df1,df2):
    common_cluster = pd.DataFrame()
    for i,j in df1.iterrows():
        list = j['genes in cluster'].split(', ')

        row =  df2[df2['genes in cluster'].str.contains('|'.join(list))]
        row1 = row.dropna()
        if row1.empty:
            pass
        else:
            data = pd.Series([row1['cluster number'].values[0],row1['genes in cluster'].values[0],row1['resource'].values[0][:-39],j['cluster number'],j['genes in cluster'],j['resource'][:-39]])
            common_cluster = common_cluster.append(data,ignore_index=True)
            outputfile = row1['resource'].values[0][:-39]+'_'+j['resource'][:-39]+'_common_clusters.tsv'
    #print(outputfile);exit(0)
    common_cluster.columns = ['cluster_number1','genes_in_cluster1','resource1','cluster_number2','genes_in_cluster2','resource2']
    common_cluster[['cluster_number1','cluster_number2']] = common_cluster[['cluster_number1','cluster_number2']].astype('int64')
    common_cluster.to_csv(outputfile,sep="\t",header=True,index=False)


common_cluster_two_dataset(read_file('Langlais_2016_gene_promoter_cluster_tf_binding.table'),read_file('Mancino_2015_gene_promoter_cluster_tf_binding.table'))
