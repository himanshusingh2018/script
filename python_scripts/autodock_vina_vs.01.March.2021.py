import pandas as pd 
import os
import pubchempy as pcp

def convert_sdf_to_pdbqt(directory):
	for ligand in os.listdir(directory):
		try:
			os.makedirs('ligand_pdbqt', exist_ok = True)
			cmd = '"c:\\Program Files (x86)\\OpenBabel-2.3.1\\babel.exe" -isdf '+directory+'\\'+ligand+' -O ligand_pdbqt\\'+ligand[:-4]+'.pdbqt'
			os.system(cmd)
			print('ligand_pdbqt\\'+ligand[:-4]+'.pdbqt is generated...')
			
		except:
			print(directory+'\\'+ligand+' is not converted into PDBQT format...')

def convert_pdb_to_pdbqt(directory):
	for ligand in os.listdir(directory):
		try:
			os.makedirs('ligand_pdbqt', exist_ok = True)
			cmd = '"c:\\Program Files (x86)\\OpenBabel-2.3.1\\babel.exe" -ipdb '+directory+'\\'+ligand+' -O ligand_pdbqt\\'+ligand[:-4]+'.pdbqt'
			os.system(cmd)
			print('ligand_pdbqt\\'+ligand[:-4]+'.pdbqt is generated...')
			
		except:
			print(directory+'\\'+ligand+' is not converted into PDBQT format...')

#convert_pdb_to_pdbqt('ligand')
			
def virtual_screening_vina(receptor_pdbqt):
	for ligand in os.listdir('ligand_pdbqt'):
		try:
			#if(TIP000001_output):
			#if not os.path.exists('./output/tf_analysis'):
			os.makedirs('docked_struct', exist_ok = True)
			if not os.path.exists('docked_struct\\'+ligand[:-6]+'_output.pdbqt'):
				cmd = '"c:\\Program Files (x86)\\PyRx\\vina.exe" --receptor '+receptor_pdbqt+' --ligand ligand_pdbqt\\'+ligand+' --config config.txt --log docked_struct\\'+ligand[:-6]+'_log.txt --out docked_struct\\'+ligand[:-6]+'_output.pdbqt'
				os.system(cmd)
				#print(cmd)
				#exit(0)
				print('Files generated:\n\tdocked_struct\\'+ligand[:-6]+'_log.txt\n\tdocked_struct\\'+ligand[:-6]+'_output.pdbqt')
			else:
				print('docked_struct\\'+ligand+'_output.pdbqt is already exist...')
				
		except:
			print('ligand_pdbqt\\'+ligand+' is not docked with receptor...')
			#exit(0)
#virtual_screening_vina('1K1B1A_BCL3.pdbqt')
	
def dock_table(docked_struct_directory):
	pdbqt = [f for f in os.listdir(docked_struct_directory) if f.endswith('_output.pdbqt')]
	df = pd.DataFrame()
	for file in pdbqt:
		with open('docked_struct\\'+file) as f:
			x=f.readlines()
			df = df.append(pd.DataFrame({'CID':[x[2].split(' = ')[1].split('\n')[0]],
								         'dG':[x[1].split('      ')[1]]}), ignore_index=True)
	#print(df.head(2))
	return(df)
dG = dock_table('docked_struct')


def drug_properties(docked_struct_directory):
	drug_properties = ['iupac_name','molecular_formula','molecular_weight',
			 		   'h_bond_acceptor_count','h_bond_donor_count','xlogp','tpsa']
	pdbqt = [f for f in os.listdir(docked_struct_directory) if f.endswith('_output.pdbqt')]
	df = pd.DataFrame()
	#os.remove('drug_properties.tsv')
	for file in pdbqt:
		with open('docked_struct\\'+file) as f:
			temp = pd.DataFrame(pcp.get_properties(drug_properties, 
									  f.readlines()[2].split(' = ')[1].split('\n')[0], 
									  'cid', as_dataframe=True))
			df = df.append(temp,ignore_index=False)
	df.reset_index(level=0, inplace=True)
	return(df)
drug_prop = drug_properties('docked_struct')

dG.to_csv('free_energy_docked_complex.tsv',sep="\t",header=True,index=False)
#print(dG)
drug_prop.to_csv('drug_properties.tsv',sep="\t",header=True,index=False)
#print(drug_prop)

df1 = pd.read_csv('free_energy_docked_complex.tsv',sep="\t",header=0)
df2 = pd.read_csv('drug_properties.tsv',sep="\t",header=0)
df = pd.merge(df1,df2, on='CID')
df.to_csv('final_complex_properties.tsv',sep="\t",header=True,index=False)
os.remove('drug_properties.tsv')
os.remove('free_energy_docked_complex.tsv')
print('final_complex_properties.tsv is generated...')
