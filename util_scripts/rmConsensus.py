#!/usr/bin/env python

'''
Remove consensus sequence (1.) from the MSAs
'''

import sys
from Bio import SeqIO

inp = sys.argv[1]
seqs = []

for idx, record in enumerate(SeqIO.parse(inp, 'fasta')):
    if idx != 0:
        seqs.append(record)

SeqIO.write(seqs, 'trunc_' + inp, 'fasta')
