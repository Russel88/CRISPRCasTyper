#!/usr/bin/env python

import multiprocessing as mp
import numpy as np
import pandas as pd
import subprocess
import re
import argparse
import os
import sys
from scipy import ndimage
from itertools import chain
from pathlib import Path
from Bio import SeqIO

# For boolean arguments
def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

ap = argparse.ArgumentParser()
ap.add_argument('-i', '--input', help='Input fasta file', required=True)
ap.add_argument('-o', '--output', help='Output prefix', required=True)
ap.add_argument('-t', '--threads', default=1, help='Number of parallel processes. Default 1', type=int)
ap.add_argument('-d', '--dist', default=1, type=int, help='Max allowed distance between genes in operon. Default 1')
ap.add_argument('-p', '--prodigal', default='single', type=str, help='Which mode to run prodigal in. Default single')
ap.add_argument('-s', '--scores', help='Path to CasScoring table', default='/mibi/users/russel/Projects/CasPredict/CasScoring.csv')
ap.add_argument('-oev', '--overall_eval', help='Overall E-value threshold. Defalt 0.001', default=0.001, type=float)
ap.add_argument('-ocs', '--overall_cov_seq', help='Overall sequence coverage threshold. Default 0.5', default=0.5, type=float)
ap.add_argument('-och', '--overall_cov_hmm', help='Overall HMM coverage threshold. Default 0.5', default=0.5, type=float)
ap.add_argument('-tev', '--two_gene_eval', help='Two-gene operon E-value threshold. Default 1e-5', default=1e-5, type=float)
ap.add_argument('-tcs', '--two_gene_cov_seq', help='Two-gene operon sequence coverage threshold. Default 0.8', default=0.8, type=float)
ap.add_argument('-tch', '--two_gene_cov_hmm', help='Two-gene operon HMM coverage threshold. Default 0.8', default=0.8, type=float)
ap.add_argument('-sev', '--single_gene_eval', help='Lonely gene E-value threshold. Default 1e-10', default=1e-10, type=float)
ap.add_argument('-scs', '--single_gene_cov_seq', help='Lonely gene sequence coverage threshold. Default 0.9', default=0.9, type=float)
ap.add_argument('-sch', '--single_cov_hmm', help='Lonely gene HMM coverage threshold. Default 0.9', default=0.9, type=float)
ap.add_argument('-vfe', '--vf_eval', help='V-F Cas12 specific E-value threshold. Default 1e-75', default=1e-75, type=float)
ap.add_argument('-vfc', '--vf_cov_hmm', help='V-F Cas12 specific HMM coverage threshold. Default 0.97', default=0.97, type=float)
ap.add_argument('-ci', '--check_input', help='Should the input be checked. Default True', default=True, type=str2bool)

# Extract arguments
args = ap.parse_args()

fasta = args.input
out = args.output
threads = args.threads
dist = args.dist
prod = args.prodigal
scoring = args.scores
oev = args.overall_eval
ocs = args.overall_cov_seq
och = args.overall_cov_hmm
tev = args.two_gene_eval
tcs = args.two_gene_cov_seq
tch = args.two_gene_cov_hmm
sev = args.single_gene_eval
scs = args.single_gene_cov_seq
sch = args.single_cov_hmm
vfe = args.vf_eval
vfc = args.vf_cov_hmm
check_inp = args.check_input

### Check input
def is_fasta(input):
    with open(input, 'r') as handle:
        fa = SeqIO.parse(handle, 'fasta')
        return any(fa)

if (not is_fasta(fasta)) and check_inp:
    sys.exit('Input file is not in fasta format')

### Check output
out = os.path.join(out, '')
try:
    os.mkdir(out)
except FileExistsError:
    sys.exit('Directory '+out+' already exists')

### Prodigal
with open(out + 'prodigal.log', 'w') as prodigal_out:
    subprocess.run(['prodigal', '-i', fasta, '-a', out+'proteins.faa', '-p', prod], stdout=prodigal_out, stderr=prodigal_out)

### Hmmer
os.mkdir(out+'hmmer')
# Define function
def hmmsearch(hmms, out):
    hmm_name = re.sub('\.hmm', '', hmms)
    with open(out+'hmmer.log', 'a') as hmmer_out:
        subprocess.run(['hmmsearch', '--domtblout', os.path.join(out+'hmmer', hmm_name+'.tab'), os.path.join('Profiles', hmms), out+'proteins.faa'], stdout=hmmer_out, stderr=hmmer_out)
# Start multiprocess
pool = mp.Pool(threads)
# Each HMM
hmm_dump = [pool.apply(hmmsearch, args=(hmms, out)) for hmms in os.listdir('Profiles')]
# Close multiprocess
pool.close()

# Extract
hmm_lst = []
for hmm_tab in os.listdir(out+'hmmer'):
    try:
        print(os.path.join(out,'hmmer',hmm_tab))
        hmm_list.append(pd.read_csv(os.path.join(out,'hmmer',hmm_tab), sep='\t', comment='#', header=None))
    except:
        pass
hmm_df = pd.concat(hmm_lst)
print(hmm_df)

sys.exit('Stop')

# Typing
typing(hmmertab, out, threads, dist, oev, ocs, och, sev, scs, sch, tev, tcs, tch, ('II-A','II-B','II-C','V-A','V-B','V-C','V-D','V-E','V-F','V-G','V-H','V-I','V-J','VI-A','VI-B1','VI-B2','VI-C','VI-D'), vfe, vfc)

