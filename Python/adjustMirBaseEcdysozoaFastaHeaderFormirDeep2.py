'''
This script the header of mature miRNA from mirBase for mirrDeep2 processing
    (ftp://mirbase.org/pub/mirbase/CURRENT/mature.fa.zip)
How to use it?
    python3 AdjustMirBaseEcdysozoaFastaHeaderFormirDeep2.py
@author Waldeyr Mendes Cordeiro da Silva
Live Long and Prosper
'''
import os
from textwrap import wrap

from Bio import SeqIO


def formatSequence(sequence, length):
    formatedSequence = ''
    listTemp = wrap(sequence, length)
    for temp in listTemp:
        formatedSequence = formatedSequence + str(temp) + '\n'
    return formatedSequence.rstrip('\n')

try:
    if os.path.exists('mirBase_mature_Ecdysozoa.fasta'):
        file = open("mirBase_mature_Ecdysozoa_for_meerDeep2.fasta", "w+")
        print("Renaming fasta headers to fit mirDeep2 constraints...")
        for record in SeqIO.parse("mirBase_mature_Ecdysozoa.fasta", "fasta"):
            header = record.description.split("|")
            newHeader = ">" + str(header[0])
            file.write(newHeader + "\n")
            file.write(formatSequence(str(record.seq), 100) + "\n")
        file.close()
    else:
        print('File not found.')
except IOError:
    print('An exception has occured.')
print("Done :-)")
