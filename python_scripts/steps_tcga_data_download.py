# STEPS TO DOWNLOAD DATA FROM TCGA DATABASE
# Website: https://portal.gdc.cancer.gov/
# Repository/
#           Data Category: transcriptome profiling
#           Data Type: Gene Expression Quantification
#           Experimental Strategy: RNA-Seq
#           Workflow Type: HTSeq-FPKM
# Cases/
#           Primary Site: Breast
#           Disese Type: ductal and lobular neoplasms (highest # of cases)
# Select the desired samples. Add to cart. Got to the Cart
#
# Add all Files to Cart/Go to Cart/Download: biospecimen, clinical,gdc_manifest
#
#           Download: Biospeciment, Clinical, Sample Sheet, Metadata, Download/Manifest
#                Manifest: (Columns)
#                       id: download folder name
#                       filename: file downloaded in the folder name 'id'
#                       md5: 128 bit md5 hashes. It is a compact digital fingerprint of a file
#                       state: validated/non-validated
#
#                Metadata: JSON file
#
#                Sample Sheet: (Columns)
#                       File ID: downloaded folder name (same in Manifest file)	
#                       File Name: downloaded file name (same in Manifest file)	
#                       Data Category: Transcriptomic profiling	
#                       Data Type: Gene Expression Quantification	
#                       Project ID: TCGA-BRCA i.e. TCGA Project	
#                       Case ID: TCGA-AC-A3QQ i.e. Case ID is One same as Sample ID
#                       Sample ID: TCGA-AC-A3QQ-01B Initial of ID is same as Case ID, just last three letters are specific samples
#                       Sample Type: Primary Tumor
#                Clinical: clinical.tsv; exposure.tsv; family_history.tsv
#                   clinical.tsv: (columns)
#                       case_id: Don't know	
#                       case_submitter_id: case id (same in Sample Sheet $ Case ID)
#
# Data Download: GDC Data Transfer Tool (https://gdc.cancer.gov/access-data/gdc-data-transfer-tool)
#
#                 Install gdc-client and run command (https://docs.gdc.cancer.gov/Data_Transfer_Tool/Users_Guide/Data_Download_and_Upload/)     
#
#		Command to download:
#				gdc-client download -m  ../gdc_manifest_6746fe840d924cf623b4634b5ec6c630bd4c06b5.txt
#           
#       Output: it download all the files from the manifest file
#
#   LOGIC TO GET FILE_NAME.CASE_ID.PRIMARY_DIAGNOSIS
#
#   Manifest File ---> (exp file name) <---             Sample Sheet
#                       Case ID/Case Submit ID <--------/
#                              ^
#                              |
#                              |
#                       Clinical ---> Primary_diagnosis
#
#   Manifest file: 'gdc_manifest_20211027_175727.txt'
#   Sample sheet: 'gdc_sample_sheet.2021-10-27.tsv'
#   clinical: 'clinical.cart.2021-10-27/clinical.tsv'

import pandas as pd
import glob
##### ANNOTATION OF SAMPLES #####
'''
#read manifest file
fManifest = pd.read_csv('gdc_manifest_20211027_204616.txt',sep="\t",header=0,usecols = ['filename'])
#read sample sheet
fSamplesheet = pd.read_csv('gdc_sample_sheet.2021-10-27.tsv',sep="\t",header=0,usecols = ['File Name','Case ID','Sample Type'])
#read clinical file
fClinical = pd.read_csv('clinical.cart.2021-10-27/clinical.tsv',sep="\t",header=0, usecols = ['case_submitter_id','primary_diagnosis'])
fClinical['primary_diagnosis'] = fClinical['primary_diagnosis'].str.split(',').str[0]
manifest_clinical = pd.merge(fSamplesheet,fClinical,left_on='Case ID',right_on='case_submitter_id')
manifest_samplesheet_clinical = pd.merge(manifest_clinical,fManifest,left_on='File Name',right_on = 'filename')
manifest_samplesheet_clinical.drop_duplicates(inplace=True)
manifest_samplesheet_clinical.to_csv('f.tsv',sep="\t",header=True,index=False)
################################

##### EXPRESSION MATRIX #####
'''
files = glob.glob('tagc_exp/*/*.FPKM.txt.gz')
#print(files)
exp = pd.read_csv(files[0],sep="\t",header=None,usecols=[0],names=['gene'])
for f in files:
    #print(f)
    name = f.split('/')[2][:-12]
    df = pd.read_csv(f,sep="\t",header=None,names=['gene',name])
    exp = pd.merge(df,exp,how='outer',on='gene')
    print('Reading file: '+f)
    
print(exp.head(4))
exp.to_csv('breast.cancer.tagc.expression.tsv',sep="\t",header=True,index=False)
    