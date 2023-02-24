'''
1. path = 'MICHELLE_0611/'
2. dirs = list of dirs
3. dir1
    all R1 in the dir merge to one file
    all R2 in the dir merge to one file
5. dir2, dir3 etc.
'''

import os, glob
idir_path = 'MICHELLE_0611/'
odir_path = 'fastq/'

for p in os.listdir(idir_path):
    r1 = glob.glob(f'{idir_path}{p}/*R1_001.fastq.gz')
    r2 = glob.glob(f'{idir_path}{p}/*R2_001.fastq.gz')

    if len(r1) == 1:
        cmd = f'cp {r1[0]} {r2[0]} {odir_path}'
        os.system(cmd)
    if len(r1)>1:
        cmd_r1 = 'cat ' + ' '.join(r1) + f' > {odir_path}' + r1[0].split('/')[-1].split('_L00')[0] + '_merge_R1_001.fastq.gz'
        cmd_r2 = 'cat ' + ' '.join(r2) + f' > {odir_path}' + r2[0].split('/')[-1].split('_L00')[0] + '_merge_R2_001.fastq.gz'
        os.system(cmd_r1)
        os.system(cmd_r2)