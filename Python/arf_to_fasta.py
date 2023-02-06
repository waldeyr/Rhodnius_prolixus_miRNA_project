'''
This script extract the header and the sequence of a mirDeep2 resulting mapped miRNA
How to use it?
    python3 arf_to_fasta.py
@author Waldeyr Mendes Cordeiro da Silva
Live Long and Prosper
'''
import os
from textwrap import wrap
import pandas as pd
from Bio import SeqIO


def formatSequence(sequence, length):
    formatedSequence = ''
    listTemp = wrap(sequence, length)
    for temp in listTemp:
        formatedSequence = formatedSequence + str(temp) + '\n'
    return formatedSequence.rstrip('\n')

arfFiles = ['../RP1G.arf', '../RP2H.arf', '.for ./RPGland.arf']

for arfFile in arfFiles:
    try:
        if os.path.exists(arfFile):
            df = pd.read_csv(arfFile, index_col=0, sep='\t', header=None)
            fastaName = str(arfFile[3:-4])+"_mapped.fasta"
            file = open(fastaName, "w+")
            for index, row in df.iterrows():
                print()
                header = f">{index}"
                sequence = row[4]
                file.write(header + "\n")
                file.write(formatSequence(sequence, 25) + "\n")
            file.close()
        else:
            print('File not found.')
    except IOError:
        print('An exception has occurred.')
print("Done :-)")
