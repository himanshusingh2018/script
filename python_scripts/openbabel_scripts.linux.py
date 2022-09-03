import pandas as pd
import os 

def split_mergeSDF_to_individualSDF(mergedSDF,sdf2dOutputDIR):
    #create folder if not exist
    os.makedirs(sdf2dOutputDIR, exist_ok = True)
    cmd = 'babel '+mergedSDF+' -O '+sdf2dOutputDIR+'/.sdf -m'
    os.system(cmd)
    print('Individual SDF files are generated...')
    
#split_mergeSDF_to_individualSDF('antiviral_updated.sdf','SDF2d')

def sdf2D_to_sdf3D(sdf2D,sdf3DoutputDIR):
    for ligand in os.listdir(sdf2D):
        try:
            os.makedirs(sdf3DoutputDIR, exist_ok = True)
            cmd = 'babel -i sdf '+sdf2D+'/'+ligand+' -O '+sdf3DoutputDIR+'/'+ligand[:-4]+'.pdbqt'
            os.system(cmd)
            print(sdf3DoutputDIR+'/'+ligand[:-4]+'.sdf is generated...')
        except:
            print(sdf3DoutputDIR+'/'+ligand+' is not converted into 3D SDF format...')
    
#sdf2D_to_sdf3D('SDF2d','SDF3d')
    
def sdf3D_to_pdbqt(SDF3d,ligPDBQToutputDIR):
	for ligand in os.listdir(SDF3d):
		try:
			os.makedirs('ligand_pdbqt', exist_ok = True)
			cmd = 'babel -isdf '+SDF3d+'/'+ligand+' -O '+ligPDBQToutputDIR+'/'+ligand[:-4]+'.pdbqt'
			os.system(cmd)
			print(ligPDBQToutputDIR+'/'+ligand[:-4]+'.pdbqt is generated...')
			
		except:
			print(ligPDBQToutputDIR+'/'+ligand+' is not converted into PDBQT format...')
            
sdf3D_to_pdbqt('SDF3d','ligand_pdbqt')
'''
def virtual_screening_vina(receptorPDBqt,ligPDBqtDIR,dockOutputDIR):
	for ligand in os.listdir(ligPDBqtDIR):
		try:
			os.makedirs(ligPDBqtDIR, exist_ok = True)
			if not os.path.exists(ligPDBqtDIR+'/'+ligand[:-6]+'_output.pdbqt'):
				cmd = '"c:\\Program Files (x86)\\PyRx\\vina.exe" --receptor '+receptor_pdbqt+' --ligand ligand_pdbqt\\'+ligand+' --config config.txt --log docked_struct\\'+ligand[:-6]+'_log.txt --out docked_struct\\'+ligand[:-6]+'_output.pdbqt'
				os.system(cmd)
                print('Files generated:\n\tdocked_struct\\'+ligand[:-6]+'_log.txt\n\tdocked_struct\\'+ligand[:-6]+'_output.pdbqt')

            else:
				print('docked_struct\\'+ligand+'_output.pdbqt is already exist...')
				
		except:
			print('ligand_pdbqt\\'+ligand+' is not docked with receptor...')
			
#virtual_screening_vina('1K1B1A_BCL3.pdbqt',)
'''