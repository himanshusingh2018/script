import pandas as pd
import numpy as np
import os
import subprocess
import datetime
import random; import statistics
from itertools import combinations

start_time = datetime.datetime.now(); print("\nThe job started at: ",start_time,"\n")

#########################################################################
#generate empty file: observed_vs_randomly_selected_genes.table
#########################################################################

class Randomization:
    def __init__(self):
        '''Read File'''

    def hg38_unique_genes_in_genome(self,coordinatefile):
        print("##### hg38 Genes coordinates from file: 'hg38.refGene_noAltChr.txt' ...")
        self.coordinatefile = coordinatefile

        refgene_col = ["bin","name","chrom","strand","txStart","txEnd","cdsStart","cdsEnd","exonCount","exonStarts","exonEnds","score","name2","cdsStartStat","cdsEndStat","exonFrames"]
        refgene_data = pd.read_csv(self.coordinatefile, sep="\t", header=None,names=refgene_col)#read hg38 refseq coordinate file without alt chr lines
        refgene_data = refgene_data[~refgene_data.chrom.str.contains('_[A-Z]') == True]#remove all alternative chromosomes from refGene

        refgene_data['txn_length'] = abs(refgene_data['txStart']-refgene_data['txEnd'])
        refgene_data_txn_length = refgene_data.groupby('name2', as_index=False).apply(pd.DataFrame.sort_values, 'txn_length', ascending = False).reset_index(drop=True)
        refgene_data['tss'] = refgene_data.groupby(['name2']).cumcount()+1
        refgene_data.loc[refgene_data['strand'] == '-',['txStart','txEnd']] = refgene_data.loc[refgene_data['strand'] == '-',['txEnd','txStart']].values

        genes_5p_3p_coord = refgene_data[['chrom','strand','txStart','txStart','txEnd','name2','tss']].copy()
        genes_5p_3p_coord1 = genes_5p_3p_coord.loc[genes_5p_3p_coord['tss'] == 1]#only tss1 is extracted

        genes_5p_coord = genes_5p_3p_coord1[['chrom','txStart','name2']].copy()
        genes_5p_coord.columns = ['chrom','txStart','txStart1','name2']

        sort_genes_5p_coord = genes_5p_coord.sort_values(['chrom','txStart'])#sort_genes_5p_coord.to_csv(cwd+'hg38.genes_5p_coord.bed',sep="\t",index=False,header=None)#print("  hg38 genes coordinate: 'hg38.genes_5p_coord.bed'\n")

        hg38_genes = list(set(refgene_data.name2))#unique gene list from hg38 geome coordinate file
        print("\t1) hg38_genes: unique genes list from hg38 RefSeq file\n\t2) sort_genes_5p_coord: tss1 (largest transcript of gene) coordinates of all genes\n")
        return hg38_genes, sort_genes_5p_coord#two variables will be input for the function: genes_random_selection

    def mm10_unique_genes_in_genome(self,coordinatefile):
        print("##### mm10 Genes coordinates from file: 'mm10.refGene_noAltChr.txt' ...")
        self.coordinatefile = coordinatefile

        refgene_col = ["bin","name","chrom","strand","txStart","txEnd","cdsStart","cdsEnd","exonCount","exonStarts","exonEnds","score","name2","cdsStartStat","cdsEndStat","exonFrames"]
        refgene_data = pd.read_csv(self.coordinatefile, sep="\t", header=None,names=refgene_col)#read hg38 refseq coordinate file without alt chr lines
        refgene_data = refgene_data[~refgene_data.chrom.str.contains('_[A-Z]') == True]#remove all alternative chromosomes from refGene

        refgene_data['txn_length'] = abs(refgene_data['txStart']-refgene_data['txEnd'])
        refgene_data_txn_length = refgene_data.groupby('name2', as_index=False).apply(pd.DataFrame.sort_values, 'txn_length', ascending = False).reset_index(drop=True)
        refgene_data['tss'] = refgene_data.groupby(['name2']).cumcount()+1
        refgene_data.loc[refgene_data['strand'] == '-',['txStart','txEnd']] = refgene_data.loc[refgene_data['strand'] == '-',['txEnd','txStart']].values

        genes_5p_3p_coord = refgene_data[['chrom','strand','txStart','txStart','txEnd','name2','tss']].copy()
        genes_5p_3p_coord1 = genes_5p_3p_coord.loc[genes_5p_3p_coord['tss'] == 1]#only tss1 is extracted

        genes_5p_coord = genes_5p_3p_coord1[['chrom','txStart','name2']].copy()
        genes_5p_coord.columns = ['chrom','txStart','txStart1','name2']

        sort_genes_5p_coord = genes_5p_coord.sort_values(['chrom','txStart'])#sort_genes_5p_coord.to_csv(cwd+'hg38.genes_5p_coord.bed',sep="\t",index=False,header=None)#print("  hg38 genes coordinate: 'hg38.genes_5p_coord.bed'\n")

        mm10_genes = list(set(refgene_data.name2))#unique gene list from hg38 geome coordinate file
        print("\t1) mm10_genes: unique genes list from hg38 RefSeq file\n\t2) sort_genes_5p_coord: tss1 (largest transcript of gene) coordinates of all genes\n")
        return mm10_genes, sort_genes_5p_coord#two variables will be input for the function: genes_random_selection

    #Lyu_HS_hg, hg38_genes, hg38_refgene_data,778, 98, 10000
    def random_selection_genes(self,resource,gene_list,genes_coordinates,induced_genes,observed_induced_genes=2,number_of_rand_selection=1):
        print('##### Choosing random genes from unique gene list in mm10/hg38 genome')
        #generate temp folder for randomly selected genes bed file
        if not os.path.exists('./randomly_selected'):
            os.mkdir('./randomly_selected')
        self.resource = resource#Article author
        self.gene_list = gene_list#list of all gene in genome
        self.genes_coordinates = genes_coordinates#coordiantes of the genes
        self.induced_genes = induced_genes#number of induced genes
        self.observed_induced_genes = observed_induced_genes#observed induced genes
        self.number_of_rand_selection = number_of_rand_selection

        randomly_selected_genes = []
        for i in range(self.number_of_rand_selection):
            selected_genes = random.sample(self.gene_list,k=self.induced_genes)
            selected_genes_coord = self.genes_coordinates[self.genes_coordinates['name2'].isin(selected_genes)]#extraction of genes tss1 coordinates on the basis of selected gene list values
            selected_genes_coord.to_csv("./randomly_selected/"+self.resource+"_randomly_selected_"+str(i+1)+"_times."+str(self.induced_genes)+"_induced_genes",sep="\t",header=False,index=False)

        print("  Files generated: \n\trandomly_selected/randomly_selected_"+str(i+1)+"_times."+str(self.induced_genes)+"_induced_genes\n")

        print("##### Induced genes distance from 5p coordinate...")

        data = selected_genes_coord.copy()
        data.columns = ['chr','5p_gene','5p_gene1','gene']

        all_distance_combination_genes = pd.DataFrame()
        for name, group in data.groupby('chr'):
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
        #print(gene_least_distant_table);exit(0)

        #all_distance_combination_genes.to_csv("./randomly_selected/randomly_selected_"+str(i+1)+"_times."+str(self.induced_genes)+"_induced_genes.all_combination_distance", sep='\t',index=False)
        gene_least_distant_table.to_csv("./randomly_selected/"+self.resource+"_randomly_selected_"+str(i+1)+"_times."+str(self.induced_genes)+"_induced_genes.least_distance", sep='\t',index=False)
        print("  Files generated: \n\trandomly_selected/"+self.resource+"_randomly_selected_"+str(i+1)+"_times."+str(self.induced_genes)+"_induced_genes.least_distance\n")

cwd = os.getcwd()+'/'
random_selection_least_distant_genes = Randomization()
hg38_genes, hg38_refgene_data = random_selection_least_distant_genes.hg38_unique_genes_in_genome('hg38.refGene.txt')
mm10_genes, mm10_refgene_data = random_selection_least_distant_genes.mm10_unique_genes_in_genome('mm10.refGene.txt')

#random_selection_of_genes_in_cluster.genes_random_selection('Viervaahra_HS',hg38_genes, hg38_refgene_data,778, 98, 1)

#arguments genes_random_selection:
#1) mm10/hg38_genes: all mm10/hg38 genes
#2) refgene_data: mm10/hg38 genes coordinates
#3) induced_genes = 65
#4) observed induced genes: 98
#5) no of random sampling: 1
#Human Data
random_selection_least_distant_genes.random_selection_genes('Biddie_Dex',mm10_genes, mm10_refgene_data,72,1)#Biddie_Dex
random_selection_least_distant_genes.random_selection_genes('Biddie_Dex_Tetracycline',mm10_genes, mm10_refgene_data,269,1)#Biddie_Dex_Tetracycline
random_selection_least_distant_genes.random_selection_genes('Brown_TNFa',hg38_genes, hg38_refgene_data,260,1)#Brown_TNFa
random_selection_least_distant_genes.random_selection_genes('Camps_Hypoxia',hg38_genes, hg38_refgene_data,212,1)#Camps_Hypoxia
random_selection_least_distant_genes.random_selection_genes('Cardamone_MitoInfec',mm10_genes, mm10_refgene_data,438,1)#Cardamone Mouse
random_selection_least_distant_genes.random_selection_genes('Demeyer_NUP214-ABL1',mm10_genes, mm10_refgene_data,5323,1)#Demeyer NUP214-ABL1
random_selection_least_distant_genes.random_selection_genes('Ebisuya_FGFstimulated',mm10_genes, mm10_refgene_data,380,1)#Ebisuya FGF-stimulated
random_selection_least_distant_genes.random_selection_genes('Ensault_Serum',mm10_genes, mm10_refgene_data,1439,1)#Ensault Serum
random_selection_least_distant_genes.random_selection_genes('Ferrari_SerumStarvation',hg38_genes, hg38_refgene_data,1274,1)#Ferrari Serum starvation
random_selection_least_distant_genes.random_selection_genes('Franco_E2',hg38_genes, hg38_refgene_data,591,1)#Franco E2 human
random_selection_least_distant_genes.random_selection_genes('Franco_TNFa',hg38_genes, hg38_refgene_data,397,1)#Franco TFNa human
random_selection_least_distant_genes.random_selection_genes('Franco_E2_TNFa',hg38_genes, hg38_refgene_data,777,1)#Franco E2+TNFa human
random_selection_least_distant_genes.random_selection_genes('Gualdrini_Starvation',mm10_genes, mm10_refgene_data,1597,1)#Gualdrini Starvation
random_selection_least_distant_genes.random_selection_genes('Hancok_Serum',hg38_genes, hg38_refgene_data,1217,1)#Hancock Serum
random_selection_least_distant_genes.random_selection_genes('Hogan_TNFa',hg38_genes, hg38_refgene_data,674,1)#Hogan TNFa human
random_selection_least_distant_genes.random_selection_genes('Hogan_IL1b',hg38_genes, hg38_refgene_data,820,1)#Hogan IL1b human
random_selection_least_distant_genes.random_selection_genes('Jin_TNFa',hg38_genes, hg38_refgene_data,613,1)#Jin TNFa human
random_selection_least_distant_genes.random_selection_genes('Jubb_GR.hMDM',hg38_genes, hg38_refgene_data,183,1)#Jubb GR.mHMDM
random_selection_least_distant_genes.random_selection_genes('Jubb_GR.mBMDM',mm10_genes, mm10_refgene_data,126,1)#Jubb GR.mBMDM
random_selection_least_distant_genes.random_selection_genes('Kusnadi_TNF',hg38_genes, hg38_refgene_data,1044,1)#Kusnadi TNF
random_selection_least_distant_genes.random_selection_genes('Langlais_IFNg',mm10_genes, mm10_refgene_data,561,1)#Langlais IFNg
random_selection_least_distant_genes.random_selection_genes('Lyu_HS',hg38_genes, hg38_refgene_data,1155,1)#Lyu HS human
random_selection_least_distant_genes.random_selection_genes('Mahat_HS',mm10_genes, mm10_refgene_data,421,1)#Mahat HS
random_selection_least_distant_genes.random_selection_genes('Mancino_LPS',mm10_genes, mm10_refgene_data,252,1)#Mancino LPS Mouse
random_selection_least_distant_genes.random_selection_genes('Park_TNF_IFN',hg38_genes, hg38_refgene_data,409,1)#Park TNF+IFN human
random_selection_least_distant_genes.random_selection_genes('Park_TNF',hg38_genes, hg38_refgene_data,997,1)#Park TNF human
random_selection_least_distant_genes.random_selection_genes('Phanstiel_PMA',hg38_genes, hg38_refgene_data,4738,1)#Phanstiel PMA
random_selection_least_distant_genes.random_selection_genes('Piccolo_IFNg',mm10_genes, mm10_refgene_data,896,1)#Piccolo IFNg Mouse
random_selection_least_distant_genes.random_selection_genes('Piccolo_IL4',mm10_genes, mm10_refgene_data,370,1)#Piccolo IL4 Mouse
random_selection_least_distant_genes.random_selection_genes('Proter_HS',hg38_genes, hg38_refgene_data,1155,1)#Porter HS human
random_selection_least_distant_genes.random_selection_genes('Schmidt_TNF',hg38_genes, hg38_refgene_data,2899,1)#Schmidt TNF
random_selection_least_distant_genes.random_selection_genes('Santiago_IFN',hg38_genes, hg38_refgene_data,482,1)#k562 encode
random_selection_least_distant_genes.random_selection_genes('Vierbuchen_Serum',mm10_genes, mm10_refgene_data,1011,1)#Vierbuchen_Serum
random_selection_least_distant_genes.random_selection_genes('Viervaahra_HS',hg38_genes, hg38_refgene_data,778,1)#Viervaahra HS human

'''
random_selection_of_genes_in_cluster.genes_random_selection('Viervaahra_HS_hg',hg38_genes, hg38_refgene_data,798,221, 1)#Viervaahra HS human
random_selection_of_genes_in_cluster.genes_random_selection('Porter_DSB_hg',hg38_genes, hg38_refgene_data,1155,138, 1)#Porter DSB human
random_selection_of_genes_in_cluster.genes_random_selection('Lyu_HS_hg',hg38_genes, hg38_refgene_data,1155,312, 1)#Lyu HS human
random_selection_of_genes_in_cluster.genes_random_selection('Jin_TNFa_hg',hg38_genes, hg38_refgene_data,613,77, 1)#Jin TNFa human
random_selection_of_genes_in_cluster.genes_random_selection('Hogan_TNFa_hg',hg38_genes, hg38_refgene_data,674, 115,1)#Hogan TNFa Human
random_selection_of_genes_in_cluster.genes_random_selection('Hogan_IL1b_hg',hg38_genes, hg38_refgene_data,820,165, 1)#Hogan IL1b human
random_selection_of_genes_in_cluster.genes_random_selection('Park_TNF_IFN_hg',hg38_genes, hg38_refgene_data,409, 78,1)#Park TNF + IFN
random_selection_of_genes_in_cluster.genes_random_selection('Park_TNF_hg',hg38_genes, hg38_refgene_data,997,230, 1)#Park TNF
random_selection_of_genes_in_cluster.genes_random_selection('Franco_E2_hg',hg38_genes, hg38_refgene_data,591,106, 1)#Franco E2 human
random_selection_of_genes_in_cluster.genes_random_selection('Franco_TNF_hg',hg38_genes, hg38_refgene_data,397, 55,1)#Franco TNFa human
random_selection_of_genes_in_cluster.genes_random_selection('Franco_E2_TNF_hg',hg38_genes, hg38_refgene_data,777, 147,1)#Franco E2+TNFa human
random_selection_of_genes_in_cluster.genes_random_selection('k562_IRF1_hg',hg38_genes, hg38_refgene_data,482, 107,1)#k562 ENCODE human
random_selection_of_genes_in_cluster.genes_random_selection('Jubb_DEX_hMDM_hg',hg38_genes, hg38_refgene_data,224, 19,1)#Jubb Dex hMDM
random_selection_of_genes_in_cluster.genes_random_selection('Kusnadi_TNF_SREBP_hg',hg38_genes, hg38_refgene_data,1115, 297,1)#Kusnadi TNF SREBP human

#Mouse Data
random_selection_of_genes_in_cluster.genes_random_selection('Mancino_LPS_mm',mm10_genes, mm10_refgene_data,252,60, 1)#Mancino LPS Mouse
random_selection_of_genes_in_cluster.genes_random_selection('Picolo_IFN_mm',mm10_genes, mm10_refgene_data,896, 273,1)#Picolo IFNg Mouse
random_selection_of_genes_in_cluster.genes_random_selection('Picolo_IL4_mm',mm10_genes, mm10_refgene_data,370,76, 1)#Picolo IL4 Mouse
random_selection_of_genes_in_cluster.genes_random_selection('Gualdrini_SRF_mm',mm10_genes, mm10_refgene_data,1693,611, 1)#Gualdrini SRF mouse
random_selection_of_genes_in_cluster.genes_random_selection('Ensault_Serum_mm',mm10_genes, mm10_refgene_data,1439,354, 1)#Ensault Serum  mouse
random_selection_of_genes_in_cluster.genes_random_selection('Jubb_Dex_mBMDM_mm',mm10_genes, mm10_refgene_data,160, 19,1)#Jubb Dex mBMDM
random_selection_of_genes_in_cluster.genes_random_selection('Cardamone_MitoInfect_mm',mm10_genes, mm10_refgene_data,438, 54,1)#Cardamone Mito Infection Mouse
'''


#os.system('rm -r ./temp')#delete all temperary files
print("The job is completed at: ",datetime.datetime.now())
print("Total time cost: ",datetime.datetime.now()-start_time)
