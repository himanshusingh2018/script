import os
with open('file.txt', 'r') as file:
    '''
    To ceate file from the paths given in file.txt
    MICHELLE_0611/Sample_14_0_1_6U_IGO_12502_L_7/14_0_1_6U_IGO_12502_L_7_S47_L003_R1_001.fastq.gz 
    '''
    for line in file:
        substring = line[:line.rfind('/')+1]#line.rfind('/'): index of last '/'
        print(substring)
        os.system(f'mkdir -p {substring}')
        os.system(f'touch {line}')

