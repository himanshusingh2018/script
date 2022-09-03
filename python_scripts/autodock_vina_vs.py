import pandas as pd 
import os

def convert_sdf_to_pdbqt(directory):
	for ligand in os.listdir(directory):
		try:
			os.makedirs('ligand_pdbqt', exist_ok = True)
			cmd = '"c:\\Program Files (x86)\\OpenBabel-2.3.1\\babel.exe" -isdf '+directory+'\\'+ligand+' -O ligand_pdbqt\\'+ligand[:-4]+'.pdbqt'
			os.system(cmd)
			print('ligand_pdbqt\\'+ligand[:-4]+'.pdbqt is generated...')
			
		except:
			print(directory+'\\'+ligand+' is not converted into PDBQT format...')
			
def virtual_screening_vina(receptor_pdbqt):
	for ligand in os.listdir('ligand_pdbqt'):
		try:
			os.makedirs('docked_struct', exist_ok = True)
			cmd = '"c:\\Program Files (x86)\\PyRx\\vina.exe" --receptor '+receptor_pdbqt+' --ligand ligand_pdbqt\\'+ligand+' --config config.txt --log docked_struct\\'+ligand[:-6]+'_log.txt --out docked_struct\\'+ligand[:-6]+'_output.pdbqt'
			os.system(cmd)
			#print(cmd)
			print('Files generated:\n\tdocked_struct\\'+ligand[:-6]+'_log.txt\n\tdocked_struct\\'+ligand[:-6]+'_output.pdbqt')
		except:
			print('ligand_pdbqt\\'+ligand+' is not docked with receptor...')
			
			

#convert_sdf_to_pdbqt('ro5')
virtual_screening_vina('prot_a_chain.pdbqt')