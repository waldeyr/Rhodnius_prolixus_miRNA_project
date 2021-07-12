import os
from Bio import SeqIO
import sys
'''
This script collapse reads in fastq format to fit the mirDeep2 input format.
It takes three arguments: the input reads, the output formated fasta file, a prefix of size=3 to name the experiment.
It has as dependencies the package Biopython version 1.78.
How touse it?
    `python3 mirDepp2_fit_fastq_to_fasta.py reads.fastq read.fasta EXP`
@author Waldeyr Mendes Cordeiro da Silva Dez-2020
Live Long and Prosper!
'''
def clear_sequence(_sequence):
    sequence = ''
    _sequence = _sequence.upper()
    for base in _sequence:
        if base in 'ATCGN':
            sequence += base
        else:
            sequence += 'N'
    return sequence

fastqFile = sys.argv[1]
fastaCollapsedFile = sys.argv[2]
readPrefix = str(sys.argv[3])
try:
    if os.path.exists(fastqFile):
        print("Parsing the fastq file...")
        records = SeqIO.parse(open(fastqFile), 'fastq')
        list_read_seq = []
        sequence_collapse = {}
        print("Collapsing reads...")
        for record in records:
            list_read_seq.append(str(record.seq).strip())
            if str(record.seq) not in sequence_collapse:
                sequence_collapse[str(record.seq)] = 0
            sequence_collapse[str(record.seq)] += 1
        print("Writing the mirDeep2 input file...")
        seqUniqueNumber = 1
        with open(fastaCollapsedFile, 'w') as writer:
            for seq, number in sequence_collapse.items():
                if len(seq) > 17:
                    # id = getMD5(seq)
                    seq = clear_sequence(seq)
                    writer.write(f">{readPrefix}_{seqUniqueNumber}_x{number}\n{seq}\n")
                    seqUniqueNumber = seqUniqueNumber + 1
        print("Done :-)")
    else:
        os.system(f'The file {fastqFile} was not found.')
except IOError:
    print('An exception occurred')