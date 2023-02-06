from Bio import SeqIO
import os
import sys


def get_ids(fname):
    reader = SeqIO.parse(fname, 'fasta')
    ids = map(lambda x: x.id, reader)
    return set(ids)

f1 = sys.argv[1]
f2 = sys.argv[2]


s1 = get_ids(f1)
s2 = get_ids(f2)

print("# Diff reads")
diff_reads = s1 - s2
for diff in diff_reads:
	print(diff)
