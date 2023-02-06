'''
This script retrieves a tab file with miRNA ID and the Species that it belongs to
    (ftp://mirbase.org/pub/mirbase/CURRENT/mature.fa.zip)
How to use it?
    python3 getMirBaseEcdysozoaSpecies.py
@author Waldeyr Mendes Cordeiro da Silva
Live Long and Prosper
'''
import os

from Bio import SeqIO

try:
    if os.path.exists('mirBase_mature_Ecdysozoa.fasta'):
        file = open("mirBase_mature_Ecdysozoa_Species.tab", "w+")
        print("Processing headers...")
        for record in SeqIO.parse("mirBase_mature_Ecdysozoa.fasta", "fasta"):
            header = record.description.split("|")
            miRNA = str(header[0])
            Species = str(header[2])
            file.write(miRNA + "\t" + Species + "\n")
        file.close()
    else:
        print('File not found.')
except IOError:
    print('An exception has occured.')
print("Done :-)")
