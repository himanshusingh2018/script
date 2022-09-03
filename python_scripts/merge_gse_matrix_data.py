import pandas as pd

def geo_exp_data_annotation(gpl_probe_annotation_file):
    '''
    Extract Probe ID and Corresponding Gene Symbol
    '''
    print('#### ---- Input File: GPL Annotation text file only ---- ####')
    #find line number where annotation of probe started
    for n,line in enumerate(open(gpl_probe_annotation_file, 'r', encoding='utf-8')):
        if "ID	GB_ACC	SPOT_ID" in line: probe_annot_row = n
    
    gpl_probe_annot = pd.read_csv(gpl_probe_annotation_file,sep="\t",skiprows=probe_annot_row,header=0,usecols = ['ID','Gene Symbol'])
    gpl_probe_annot['Gene Symbol'] = gpl_probe_annot['Gene Symbol'].str.split(' ').str[0]
    print(gpl_probe_annot.head(4))
    return(gpl_probe_annot)


def read_geo_matrix(geo_matrix_file,probe_annot):
    '''
    Extract all expression data from the Microarray matrix dataset of specific GSE ID
    '''
    print('#### ---- Input File: GEO Matrix text file only ---- ####')
    #find line number where annotation of samples given and where the data started
    for n,line in enumerate(open(geo_matrix_file, 'r', encoding='utf-8')):
        if "ID_REF" in line: nskiprows = n
        if "!Sample_title" in line: sample_annot_row = n-1
    
    sample_annot = pd.read_csv(geo_matrix_file, sep="\t",header=0,skiprows = sample_annot_row,nrows=0)
    exp = pd.read_csv(geo_matrix_file,sep="\t",header=None,skiprows=nskiprows,names=sample_annot.columns)
    probe_annot = probe_annot.append(pd.DataFrame({'ID':['ID_REF'],'Gene Symbol':['GSM_ID']})).reset_index(drop=True)
    probe_annot.rename(columns={'ID':'!Sample_title'}, inplace=True)
    annot_exp = pd.merge(probe_annot, exp, how="outer", on=["!Sample_title"])
    annot_exp.to_csv(geo_matrix_file.split('_')[0]+'_expression.txt',sep="\t",header=True,index=False)
    print('Expression file successfully generated:\n\t'+geo_matrix_file.split('_')[0]+'_expression.txt')

read_geo_matrix(geo_matrix_file = 'GSE29532_series_matrix.txt', probe_annot = geo_exp_data_annotation(gpl_probe_annotation_file = 'GPL5175-3188.txt'))
#read_geo_matrix(geo_matrix_file = 'GSE62646_series_matrix.txt', probe_annot = geo_exp_data_annotation(gpl_probe_annotation_file = 'GPL6244-17930.txt'))
#read_geo_matrix(geo_matrix_file = 'GSE66360_series_matrix.txt', probe_annot = geo_exp_data_annotation(gpl_probe_annotation_file = 'GPL570-55999.txt'))

    
