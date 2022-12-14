import datetime
import pandas as pd
import numpy as np
import os
import subprocess
import seaborn as sns
from itertools import combinations
import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager

os.chdir(os.getcwd())#moving to current working directory
start_time = datetime.datetime.now()
print("\nThe job started at: ",datetime.datetime.now())

input_files ="""
All the required files must be avilable in folder: 'input_data'\n
Required Input files:
1: Expression Data file required: 'Induced.genes'
   file with at least two below header line
   [gene  log2FC]\n
2: RefSeq GTF file: 'hg38.refGene_noAltChr.txt'
   [bin,name,chrom,strand,txStart,txEnd,cdsStart,cdsEnd,exonCount,exonStarts,exonEnds,score,name2,cdsStartStat,cdsEndStat,exonFrames]\n
3: PeakFile: peak.narrowPeak
   with columns without header line. No header line in peak file. Below columns must be there
   ['chr','start','end','peak.name','peak.size','strand']\n
4: Filter Score: default: 0, User can define from: peak.narrowPeak
   ['chr', 'start', 'end', 'length', 'abs_summit', 'pileup', '-log10(pvalue)', 'fold_enrichment', '-log10(qValue)', 'name']
"""
print(input_files)

if not os.path.exists('./output/tf_analysis'):#generate output folders: exp_bedgraph and tf_analysis
    for sub_folders in ['output/exp_bedgraph','output/tf_analysis']:
        os.makedirs(sub_folders)

class Expression_bedgraph:
    '''Expression_bedgraph'''
    def __init__(self):
        '''Read File'''
    def induced_genes(self,expression_file,gtf):
        os.chdir(cwd)
        print("##### Generation of bedgraph from the induced genes...")
        self.expression_file = expression_file
        self.gtf = gtf
        #expression = pd.read_excel('input_data/'+self.expression_file,skiprows=7,sheet_name='table S1B',header=0,sep="\t",usecols=['Official Gene Symbol','logFC'])
        expression = pd.read_csv('input_data/'+self.expression_file,header=0,sep='\t')
        expression['log2FC'] = (expression.logCPM_T1 + expression.logCPM_T2 + expression.logCPM_T3)-(expression.logCPM_V1 + expression.logCPM_V2 + expression.logCPM_V3)

        expression = expression[['Symbol','log2FC']]
        expression.columns = ['name2','DESeq2_log2FoldChange']

        induced_gene_exp_data = expression[(expression.DESeq2_log2FoldChange > 1)]
        #writing induced_genes
        induced_gene_exp_data.to_csv('./output/exp_bedgraph/'+self.expression_file.split('.')[0]+'.hg19.induced_genes',sep="\t",header=False,index=False)

        refgene_col = ["bin","name","chrom","strand","txStart","txEnd","cdsStart","cdsEnd","exonCount","exonStarts","exonEnds","score","name2","cdsStartStat","cdsEndStat","exonFrames"]
        refgene_data = pd.read_csv('input_data/'+self.gtf,sep="\t",header=None,names=refgene_col)
        refgene_data = refgene_data[~refgene_data.chrom.str.contains('_')]
        refgene_data['txn_length'] = abs(refgene_data['txStart']-refgene_data['txEnd'])
        refgene_data_txn_length = refgene_data.groupby('name2', as_index=False).apply(pd.DataFrame.sort_values, 'txn_length', ascending = False).reset_index(drop=True)
        induced_refgene_data = refgene_data_txn_length[refgene_data_txn_length['name2'].isin(induced_gene_exp_data['name2'])].copy()
        induced_refgene_data['tss'] = induced_refgene_data.groupby(['name2']).cumcount()+1

        exp_with_coordinate = induced_gene_exp_data.merge(induced_refgene_data, on='name2')
        exp_with_coordinate_all_tss = exp_with_coordinate[["name2", "DESeq2_log2FoldChange", "chrom","strand", "txStart","txEnd","tss"]].copy()
        #exp_with_coordinate_tss1 = exp_with_coordinate_all_tss.loc[exp_with_coordinate_all_tss['tss'] == 1]

        exp_bedgraph = pd.DataFrame()
        for item,row in exp_with_coordinate_all_tss.iterrows():
            if(row['strand'] == '+'):
                exp_bedgraph = exp_bedgraph.append(pd.Series([row['chrom'],row['txStart'],row['txStart']+500,round(row['DESeq2_log2FoldChange'],4),row['name2']]),ignore_index=True)
            if(row['strand'] == '-'):
                exp_bedgraph = exp_bedgraph.append(pd.Series([row['chrom'],row['txEnd'],row['txEnd']-500,round(row['DESeq2_log2FoldChange'],4),row['name2']]),ignore_index=True)

        exp_bedgraph.columns = ['chr','sTSS','TSSto500bps','log2FC','geneName']
        exp_bedgraph = exp_bedgraph.astype({"sTSS": int, "TSSto500bps": int})
        exp_bedgraph.to_csv('output/exp_bedgraph/'+self.expression_file.split('.')[0]+'.hg19.bedgraph',sep="\t",header=False,index=False)
        print(' Expression Bedgraph file is generated:\n\toutput/exp_bedgraph/'+self.expression_file.split('.')[0]+'.hg19.bedgraph',"\n")

class Extract_induced_genes_coordinates:
    def __init__(self):
        '''Read File'''
    def induced_genes(self,genelist,coordinatefile):
        print("##### Extraction induced genes coordinates from RefSeq Database file: 'RefGene.txt' ...")
        self.coordinatefile = coordinatefile
        self.genelist = genelist

        induced_gene_data = pd.read_csv('output/exp_bedgraph/'+self.genelist[:-4]+'.hg19.induced_genes',sep="\t",header=0)
        induced_gene_data.columns = ['gene','log2FC']

        refgene_col = ["bin","name","chrom","strand","txStart","txEnd","cdsStart","cdsEnd","exonCount","exonStarts","exonEnds","score","name2","cdsStartStat","cdsEndStat","exonFrames"]
        refgene_data = pd.read_csv('input_data/'+self.coordinatefile, sep="\t", header=None,names=refgene_col)
        refgene_data['txn_length'] = abs(refgene_data['txStart']-refgene_data['txEnd'])
        refgene_data_txn_length = refgene_data.groupby('name2', as_index=False).apply(pd.DataFrame.sort_values, 'txn_length', ascending = False).reset_index(drop=True)

        induced_refgene_data = refgene_data_txn_length[refgene_data_txn_length['name2'].isin(induced_gene_data['gene'])].copy()
        induced_refgene_data['tss'] = induced_refgene_data.groupby(['name2']).cumcount()+1
        induced_gene_tss = pd.DataFrame(induced_refgene_data, columns=['name','chrom','strand','txStart','txEnd','name2','tss'])

        print("  Induced genes tss coordinates are saved in variable: 'induced_gene_tss'\n")
        induced_gene_tss.to_csv('output/tf_analysis/induced_gene_tss',sep="\t",index=False)
        #print(induced_gene_tss.head(4))
        return induced_gene_tss

class Inducedgenesclustering:
    def __init__(self):
        '''Clustering on the basis of 5' coordinate of induced genes'''
    def induced_genes_clustering(self, data):
        print("##### Induced genes clustering...")
        self.data = data
        os.chdir(cwd)
        five_prime_coordinate_gene = pd.DataFrame()

        dtypes = {'chr':'str','strand':'str','5p_gene1':'int','5p_gene2':'int','3p_gene':'int','gene':'str'}
        five_prime_coordinate_gene = pd.DataFrame(columns=['chr','5p_gene1','5p_gene2','3p_gene','gene'])#gene coordinate

        for index, row in self.data.iterrows():
            if(row['strand'] == '+' and row['tss'] == 1):
                five_prime_coordinate_gene=five_prime_coordinate_gene.append(pd.Series([row['chrom'],row['txStart'],
                row['txStart'],row['txEnd'],row['name2']],index=['chr','5p_gene1','5p_gene2','3p_gene','gene']),ignore_index=True)

            elif(row['strand'] == '-' and row['tss']==1):
                five_prime_coordinate_gene=five_prime_coordinate_gene.append(pd.Series([row['chrom'],row['txEnd'],row['txEnd'],
                row['txStart'],row['name2']],index=['chr','5p_gene1','5p_gene2','3p_gene','gene']),ignore_index=True)

        sort_five_prime_coordinate_gene = five_prime_coordinate_gene.sort_values(['chr','5p_gene1','5p_gene2'])
        five_prime_coordinate = sort_five_prime_coordinate_gene[['chr','5p_gene1','5p_gene2','gene']]
        five_prime_coordinate.to_csv('output/tf_analysis/induced_genes_5p_coordinates.bed',sep="\t",index=False,header=None)
        sort_five_prime_coordinate_gene.to_csv('output/tf_analysis/induced_genes_5p_3p_coordinates.bed',sep='\t',index=False,header=None)

        cmd_cluster_5p_genes = "bedtools cluster -d 100000 -i output/tf_analysis/induced_genes_5p_3p_coordinates.bed > output/tf_analysis/induced_genes_bed_tool_clust.bed"
        subprocess.call(cmd_cluster_5p_genes,shell=True)
        #count cluster which has more than one genes
        cmd_count_cluster= "awk 'NR == FNR {CNT[$NF]++;if (!($NF in order) && CNT[$NF]>1) order[$NF]=++cnt;next} $NF in order {print $0, order[$NF]}' OFS='\t' output/tf_analysis/induced_genes_bed_tool_clust.bed output/tf_analysis/induced_genes_bed_tool_clust.bed > output/tf_analysis/induced_genes_final_cluster.bed"
        subprocess.call(cmd_count_cluster, shell=True)

        print("  Induced Gene Cluster file is generated:\n\t output/tf_analysis/induced_genes_final_cluster.bed\n")
        return sort_five_prime_coordinate_gene

class Genesclusterbarplot:
    def __init__(self):
        '''Save Plot'''
    def genes_clust_genes_frequency(self, filename):
        print("##### Genes cluster bar plot...")
        self.filename = filename
        cluster_col_names = ['chr','5p_gene1','5p_gene2','3p_gene','gene','clust_no','cluster']

        hist_table_data = pd.read_csv('output/tf_analysis/'+self.filename, sep="\t", header=None,names=cluster_col_names)
        hist_table = pd.DataFrame(hist_table_data.groupby(['cluster']).size().value_counts().reset_index(name='number of clusters'))

        ax=sns.barplot(x="index",y="number of clusters",data=hist_table,color='blue')
        ax.set_xlabel("genes in cluster",fontsize=14)   #ax.axes.set_title("Title",fontsize=16)
        ax.set_ylabel("frequency of clusters",fontsize=14)
        plt.savefig("output/tf_analysis/cluster_genes_frequency_at_100kb_plot.png")

        hist_table.to_csv('output/tf_analysis/cluster_gene_frequency_at_100kb_data.txt',sep="\t",header=['genes in cluster','frequency of clusters'],index=False)
        print("  Frequency of Induced Genes per Cluster data file is generated:\n\toutput/tf_analysis/cluster_gene_frequency_at_100kb_data.txt")
        print("  Frequency of Induced Genes per Cluster Barplot is generated:\n\toutput/tf_analysis/cluster_gene_frequency_at_100kb_barplot.png\n")

class Genesleastdistant:
    def __init__(self):
        '''Save Plot'''

    def genes_least_distnce(self, filename):
        print("##### Induced genes distance from 5p coordinate...")
        self.filename = filename
        data = pd.read_csv('output/tf_analysis/'+self.filename, sep="\t",header=None,names=['chr','5p_gene','5p_gene1','3p_gene','gene'])

        all_distance_combination_genes = pd.DataFrame()
        for name, group in data.groupby('chr'):
            #print(name)
            gene_pairs = list(combinations(group.gene,2))
            dist_pairs = list(combinations(group['5p_gene'],2))
            gene_combination = pd.DataFrame(gene_pairs,columns=['gene1','gene2'])
            dist_combination = pd.DataFrame(dist_pairs,columns=['gene1_dist','gene2_dist'])

            df = pd.concat([gene_combination,dist_combination],axis=1)
            df['distance'] = abs(df.gene1_dist-df.gene2_dist)
            all_distance_combination_genes = all_distance_combination_genes.append(df[['gene1','gene2','distance']],ignore_index=True,sort=False)

        all_distance_combination_genes['distance'] = all_distance_combination_genes['distance'].astype('int64')
        gene_least_distant_genes_indices = all_distance_combination_genes.groupby('gene1')['distance'].idxmin
        gene_least_distant_table = all_distance_combination_genes.loc[gene_least_distant_genes_indices]

        title_font = {'fontname':'sans-serif', 'size':'16', 'color':'black', 'weight':'normal','verticalalignment':'bottom'} # Bottom vertical align
        axis_font = {'fontname':'sans-serif', 'size':'14'}
        bar_plot, ax = plt.subplots(nrows=1,ncols=1)
        x_obj = ["???100 kb","101-200 kb", "201-500 kb", ">500 kb"]
        y_obj = [sum(i < 100001 for i in gene_least_distant_table['distance']),
                sum(100000 < i < 200001 for i in gene_least_distant_table['distance']),
                sum(200000 < i < 500001 for i in gene_least_distant_table['distance']),
                sum(i > 500000 for i in gene_least_distant_table['distance'])]

        plt.bar(x_obj,y_obj, align = 'center')#BAR PLOT
        #plt.plot(x_obj,y_obj)#, align = 'center') #LINE PLOT
        plt.xlabel('Least distance between two genes from 5\' coordinate', **axis_font)
        plt.ylabel('Frequency of genes', **axis_font)
        bar_plot.savefig('output/tf_analysis/gene_least_distant_barplot.png')
        plt.close(bar_plot)

        all_distance_combination_genes.to_csv('output/tf_analysis/dist_all_induced_genes.distance', sep='\t',index=False)
        gene_least_distant_table.to_csv('output/tf_analysis/least_distant_induced_genes.distance', sep='\t',index=False)
        print("  Files generated: \n\toutput/tf_analysis/least_distant_induced_genes.distance\n\toutput/tf_analysis/dist_all_induced_genes.distance\n\toutput/tf_analysis/gene_least_distant_barplot.png\n")

class Tf_peak_tss_binding:
    def __init__(self):
        '''TF peak tss binding'''
    def tf_peak_tss_cluster(self,filename,tf,object,filterscore):
        os.chdir(cwd)
        print("#####TF Peak Tss Binding...")
        self.filename = filename
        self.tf = tf
        self.object = object
        self.filterscore = filterscore

        #tf_peak_cols = ['chr','start','end','peak.size','score','strand','fold_change','-log10.pValue','-log10.qValue','summit']
        #tf_peak_cols=['chr','start','end','peak.name','peak.size','strand']
        tf_peak_file = pd.read_csv('input_data/'+self.filename,sep="\t",header=0,usecols=['Chr','Start','End','Peak Score'])        #macs2 $1,$2,$3,$5
        tf_peak_file = tf_peak_file[['Chr','Start','End','Peak Score']]
        tf_peak_file.columns = ['chr','start','end','score']        #macs2 $1,$2,$3,$5

        tf_peak_file['score'] = tf_peak_file['score'].apply(pd.to_numeric, downcast='float', errors='coerce')
        filter = tf_peak_file['score'] > self.filterscore

        tf_peak_data = pd.DataFrame(tf_peak_file[filter], columns=['chr','start','end', 'score'])
        tf_peak_data['tf'] = self.tf
        tf_peak_data.loc[tf_peak_data['end'] < tf_peak_data['start'],['start','end']] = tf_peak_data.loc[tf_peak_data['end'] < tf_peak_data['start'],['end','start']].values
        gene_tss=self.object[['chrom','strand','txStart','txEnd','name2','tss']].copy()
        #swap column of negative strand
        strand = gene_tss['strand'] == '-'; gene_tss.loc[strand, ['txStart','txEnd']] = gene_tss.loc[strand,['txEnd','txStart']].values
        del gene_tss['strand']#delete the strand column from dataframe
        gene_tss_5p = gene_tss[['chrom','txStart','txStart','name2','tss']].copy()
        #print(gene_tss_5p.head(2),"\n\n",gene_tss.head(2));exit(0)
        gene_tss_5p.columns = tf_peak_data.columns

        merge_tf_tss = pd.concat([tf_peak_data,gene_tss_5p])
        merge_tf_gene = merge_tf_tss.sort_values(['chr','start','end'])

        merge_tf_gene.loc[merge_tf_gene['end'] < merge_tf_gene['start'],['start','end']] = merge_tf_gene.loc[merge_tf_gene['end'] < merge_tf_gene['start'],['end','start']].values
        sort_merge_tf_gene = merge_tf_gene.sort_values(['chr','start'])
        #print(sort_merge_tf_gene.head(2));exit(0)

        sort_merge_tf_gene['start'] = sort_merge_tf_gene['start'].astype(int)
        sort_merge_tf_gene['end'] = sort_merge_tf_gene['end'].astype(int)
        #print(sort_merge_tf_gene.head(2));exit(0)
        sort_merge_tf_gene.to_csv('output/tf_analysis/combined_gene_5p_tf_peak_calls_data.bed',sep="\t",index=False,header=None)

        cmd_tss_peak_clust = "bedtools cluster -d 1000 -i output/tf_analysis/combined_gene_5p_tf_peak_calls_data.bed > output/tf_analysis/gene_5p_tf_peak_calls_bed_tool_clust.bed"
        subprocess.call(cmd_tss_peak_clust,shell=True)
        print("  Gene Tss and TF Peak Cluster file is generated:\n\t"+cwd+"gene_5p_tf_peak_calls_bed_tool_clust.bed")

        cluster_file_columns = ['chr', 'sTSS', 'eTSS', 'gene', 'tf', 'clust']
        cluster_file = pd.read_csv('output/tf_analysis/gene_5p_tf_peak_calls_bed_tool_clust.bed', sep='\t', names=cluster_file_columns, dtype={'clust_no':int, 'tf_gene_tss_distance': int}, low_memory=False)

        tf_df = cluster_file[cluster_file['tf'] == self.tf].reset_index(drop=True)
        non_tf_df = cluster_file[cluster_file['tf'] != self.tf].reset_index(drop=True)

        #print(non_tf_df.head(4));exit(0)
        temp_peak_table = pd.merge(non_tf_df, tf_df, on='clust')

        temp_peak_table.columns=['chr_gene','sTSS_gene','eTSS_gene','gene','tss_no','clust_no','chr_tf','sTSS_tf','eTSS_tf','score','tf']

        peak_table=temp_peak_table[['clust_no','gene','tss_no','score']].copy()
        peak_table['distance'] = np.minimum(abs(temp_peak_table['sTSS_gene']-temp_peak_table['eTSS_tf']),abs(temp_peak_table['sTSS_gene']-temp_peak_table['sTSS_tf']))
        #print(peak_table);exit(0)
        peak_table.to_csv("output/tf_analysis/gene_tss_peak.table", sep='\t',index=False)
        print("  The gene tss peak table is generated: \n\toutput/tf_analysis/gene_tss_peak.table\n")

class Genespromotercluster:
    def __init__(self):
        '''Gene Promoter Cluster Table'''

    def genes_promoter_cluster(self,filename):
        print("#####Genes promoter cluster...")
        self.filename = filename

        columns = ['chr', 'gene_5p', 'gene_5p2', 'gene_3p', 'gene_name', 'clust1', 'clust2']
        gene_cluster_data = pd.read_csv('output/tf_analysis/'+self.filename, sep="\t", names=columns)

        cluster_number = pd.Series(gene_cluster_data['clust2'].unique())
        gene_promoter_cluster_table_columns = ['cluster number', 'number of genes', 'number of promoters', 'flag','genes of cluster']
        gene_promoter_cluster_table = pd.DataFrame()

        for i in cluster_number:
            number_of_genes = []
            cluster_genes = []
            for j in range(0, len(gene_cluster_data)-1):
                if(gene_cluster_data.iloc[j, 6] == i):
                    k = j + 1
                    if(gene_cluster_data.iloc[k,6] == i):
                        dist = abs(gene_cluster_data.iloc[j, 1]-gene_cluster_data.iloc[k, 1])
                        number_of_genes.append(dist)
                    cluster_genes.append(gene_cluster_data.iloc[k-1, 4])
                    if(k == len(gene_cluster_data)-1):
                        cluster_genes.append(gene_cluster_data.iloc[k, 4])
            number_of_genes_per_cluster = len(number_of_genes)+1
            number_of_promoters = sum(i > 1000 for i in number_of_genes)+1
            flag = 'U' if number_of_promoters > 1 else 'F'
            gene_promoter_cluster_table=gene_promoter_cluster_table.append(pd.Series([i,number_of_genes_per_cluster, number_of_promoters, flag, ', '.join(cluster_genes)], index=gene_promoter_cluster_table_columns), ignore_index=True)

        gene_promoter_cluster_table['cluster number'] = gene_promoter_cluster_table['cluster number'].astype(int)
        gene_promoter_cluster_table['number of genes'] = gene_promoter_cluster_table['number of genes'].astype(int)
        gene_promoter_cluster_table['number of promoters'] = gene_promoter_cluster_table['number of promoters'].astype(int)
        gene_promoter_cluster_table['flag'] = gene_promoter_cluster_table['flag'].astype(str)
        gene_promoter_cluster_table['genes of cluster'] = gene_promoter_cluster_table['genes of cluster'].astype(str)
        gene_promoter_cluster_table = gene_promoter_cluster_table.reindex(columns=gene_promoter_cluster_table_columns)

        gene_promoter_cluster_table.to_csv('output/tf_analysis/gene_promoter_cluster.table', sep='\t',index=False)
        print("  The Gene Promoter Cluster Table is generated:\n\toutput/tf_analysis/gene_promoter_cluster.table\n")
        return gene_promoter_cluster_table

class Genespromotertfbindingcluster:
    def __init__(self):
        '''Gene Promoter TF Binding'''

    def tss_tf_binding_cluster(self,filename,object):
        print("#####Gene Promoter TF binding cluster...")
        self.filename = filename
        self.object = object

        gene_promoter_cluster_data = self.object.copy()
        gene_peak_data = pd.read_csv('output/tf_analysis/'+self.filename,sep="\t",header=0)
        tf_bound_genes = gene_peak_data['gene'].unique()

        gene_promoter_cluster_tf_binding_table_columns = ['cluster number','number of genes per cluster','number of promoters per cluster','flag','genes in cluster','number of promoter binding TFs in cluster']
        gene_promoter_cluster_tf_binding_table = pd.DataFrame(columns=gene_promoter_cluster_tf_binding_table_columns)


        for index, gene_promoter_cluster_data_row in gene_promoter_cluster_data.iterrows():
            count = 0

            genes_in_cluster = [x.strip() for x in gene_promoter_cluster_data_row[4].split(',')]

            for element in genes_in_cluster:
                if element in tf_bound_genes:
                    count += 1
            gene_promoter_cluster_tf_binding_table = gene_promoter_cluster_tf_binding_table.append(pd.Series([gene_promoter_cluster_data_row[0],gene_promoter_cluster_data_row[1],gene_promoter_cluster_data_row[2],gene_promoter_cluster_data_row[3],gene_promoter_cluster_data_row[4],count],index=gene_promoter_cluster_tf_binding_table_columns), ignore_index = True)

        gene_promoter_cluster_tf_binding_table.to_csv('output/tf_analysis/gene_promoter_cluster_tf_binding.table', sep='\t',index=False)

        gene_promoter_cluster_tf_binding_frequency_table = gene_promoter_cluster_tf_binding_table[gene_promoter_cluster_tf_binding_table.apply(lambda x: x['flag'] == 'U', axis=1)].groupby(['number of promoters per cluster', 'number of promoter binding TFs in cluster']).size().reset_index(name="Frequency")

        bubble_chart_frequency_table = gene_promoter_cluster_tf_binding_frequency_table.groupby(['number of promoters per cluster', 'number of promoter binding TFs in cluster'])['Frequency'].sum().reset_index()
        promoter_per_cluster = bubble_chart_frequency_table.groupby(['number of promoters per cluster'])['Frequency'].sum().reset_index()
        bubble_chart_frequency_table = bubble_chart_frequency_table.merge(promoter_per_cluster, left_on='number of promoters per cluster', right_on='number of promoters per cluster', how='left')
        bubble_chart_frequency_table['Frequency Percentage'] = 100*(bubble_chart_frequency_table['Frequency_x']/bubble_chart_frequency_table['Frequency_y'])
        bubble_chart_frequency_table.to_csv('output/tf_analysis/bubble_chart_gene_promoter_cluster_tf_binding_frequency.table',sep="\t",index=None)

        print("  The files generated:\n\toutput/tf_analysis/gene_promoter_cluster_tf_binding.table\n\toutput/tf_analysis/bubble_chart_gene_promoter_cluster_tf_binding_frequency.table\n")

class TfBindingGenePromoterFrequencyBubblePlot:
    def __init__(self):
        '''Bubble Plot'''

    def bubble_plot(self):
        print("##### TF Binding Gene Promoter Frequency Bubble Plot Generating...")
        #self.filename = filename
        cmd = '/usr/bin/Rscript tf_binding_gene_promoter_frequency_bubble_plot.10.04.2020.R'
        os.system(cmd)
        print("  The files generated:\n\toutput/tf_analysis/TF_binding_gene_promoter_frequency_bubble_plot.png\n")

cwd = os.getcwd()+'/'

#FILE NAMES PRESENT IN INPUT_DATA DIRECTORY MUST BE DECIDED
#mouse_genome_annotation_refseq_gtf = 'mm9.refGene.txt'
hg19_refseq_gtf = 'hg19.refGene.txt'

rna_expression_file = 'GSE60462_RNA_Exon_counts.txt'
tf_name = 'p65'

tf_chiqseq_bedfile1 = 'GSE64233_MED1.peakfile.txt'
tf_chiqseq_bedfile2 = 'GSE64233_MED1_peakfile.txt'
tf_chiqseq_bedfile3 = 'GSE64233_p65.peakfile.txt'


'''
TAKING INPUT FROM USER
genome_annotation_refseq_gtf = input('Enter the name of genome annotation file:')
rna_expression_file = input('Enter the name of RNA expression file:')
tf_name = input('Enter the name of TF:')
tf_chiqseq_bedfile = input('Enter the name of TF chipseq bed file:')
'''

expression_bedgraph=Expression_bedgraph()
#expression_bedgraph.induced_genes(rna_expression_file,hg19_refseq_gtf)

extract_induced_genes=Extract_induced_genes_coordinates()
induced_genes_tss = extract_induced_genes.induced_genes(rna_expression_file,hg19_refseq_gtf)

induced_genes_clust = Inducedgenesclustering()
induced_genes_cluster=induced_genes_clust.induced_genes_clustering(induced_genes_tss)

induced_genes_cluster_barplot=Genesclusterbarplot()
induced_genes_cluster_barplot.genes_clust_genes_frequency('induced_genes_final_cluster.bed')

least_distant_induced_genes= Genesleastdistant()
least_distant_induced_genes.genes_least_distnce('induced_genes_5p_3p_coordinates.bed')

tf_peak_tss_binding = Tf_peak_tss_binding()
tf_peak_tss_binding.tf_peak_tss_cluster(tf_chiqseq_bedfile3,tf_name,induced_genes_tss,filterscore=0)

genes_promoter_cluster = Genespromotercluster()
genes_promoter_clust = genes_promoter_cluster.genes_promoter_cluster('induced_genes_final_cluster.bed')

tss_tf_clustering = Genespromotertfbindingcluster()
tss_tf_cluster = tss_tf_clustering.tss_tf_binding_cluster('gene_tss_peak.table',genes_promoter_clust)

freq_bubble_plot = TfBindingGenePromoterFrequencyBubblePlot()
freq_bubble_plot.bubble_plot()


print("The job is finished by: ",datetime.datetime.now())
print("Total time cost: ",datetime.datetime.now()-start_time,"\n")
