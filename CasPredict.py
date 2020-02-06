#!/usr/bin/env python

import numpy as np
import pandas as pd
import subprocess
import glob
import re
import argparse
import logging
import os
import sys
import shutil
from pathlib import Path
from Bio import SeqIO
import multiprocessing as mp

from cas_typing import typing, type_operon, cluster_adj

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

# Required
ap.add_argument('input', help='Input fasta file')
ap.add_argument('output', help='Prefix for output directory')

# Optional
ap.add_argument('-t', '--threads', default=4, help='Number of parallel processes. Default 4', type=int)
ap.add_argument('--prodigal', default='single', type=str, help='Which mode to run prodigal in. Default single')
ap.add_argument('--aa', help='Input is a protein fasta. Has to be in prodigal format', action='store_true')
ap.add_argument('--check_input', help='Should the input be checked. Default True', default=True, type=str2bool)
ap.add_argument('--keep_prodigal', help='Keep prodigal output. Default False', default=False, type=str2bool)
ap.add_argument('--log_lvl', help='Logging level. Default 20', default=20, type=int)
ap.add_argument('--redo_typing', help='Redo the typing. Skip prodigal and HMMER and load the hmmer.tab from the output dir', action='store_true')

# Data
apd = ap.add_argument_group('data arguments')
apd.add_argument('--scores', help='Path to CasScoring table. Default same dir as CasPredict script', default='', type=str)
apd.add_argument('--hmms', help='Path to directory with HMM profiles. Default same dir as CasPredict script', default='', type=str)

# Thresholds
apt = ap.add_argument_group('threshold arguments')
apt.add_argument('--dist', default=3, type=int, help='Max allowed distance between genes in operon. Default 3')
apt.add_argument('--overall_eval', help='Overall E-value threshold. Defalt 0.001', default=0.001, type=float)
apt.add_argument('--overall_cov_seq', help='Overall sequence coverage threshold. Default 0.5', default=0.5, type=float)
apt.add_argument('--overall_cov_hmm', help='Overall HMM coverage threshold. Default 0.5', default=0.5, type=float)
apt.add_argument('--two_gene_eval', help='Two-gene operon E-value threshold. Default 1e-5', default=1e-5, type=float)
apt.add_argument('--two_gene_cov_seq', help='Two-gene operon sequence coverage threshold. Default 0.8', default=0.8, type=float)
apt.add_argument('--two_gene_cov_hmm', help='Two-gene operon HMM coverage threshold. Default 0.8', default=0.8, type=float)
apt.add_argument('--single_gene_eval', help='Lonely gene E-value threshold. Default 1e-10', default=1e-10, type=float)
apt.add_argument('--single_gene_cov_seq', help='Lonely gene sequence coverage threshold. Default 0.9', default=0.9, type=float)
apt.add_argument('--single_cov_hmm', help='Lonely gene HMM coverage threshold. Default 0.9', default=0.9, type=float)
apt.add_argument('--vf_eval', help='V-F Cas12 specific E-value threshold. Default 1e-75', default=1e-75, type=float)
apt.add_argument('--vf_cov_hmm', help='V-F Cas12 specific HMM coverage threshold. Default 0.97', default=0.97, type=float)

# Extract arguments
args = ap.parse_args()

fasta = args.input
out = args.output
threads = args.threads
dist = args.dist
prod = args.prodigal
scoring = args.scores
profile_dir=args.hmms
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
keep_prodigal = args.keep_prodigal
lvl = args.log_lvl
aa = args.aa
redo = args.redo_typing

# Force argument consistency with protein input
if aa:
    keep_prodigal = True

# Force arguments with redo
if redo:
    aa = True
    check_inp = False

# Logger
logging.basicConfig(format='\033[36m'+'[%(asctime)s] %(levelname)s:'+'\033[0m'+' %(message)s', datefmt='%Y-%m-%d %H:%M:%S', level=lvl)

# Version
logging.info('Running CasPredict version 0.1.3')

# Data dir
script_dir = re.sub('CasPredict.py', '', os.path.realpath(__file__))
if scoring == '':
    scoring = os.path.join(script_dir, 'CasScoring.csv')
if profile_dir == '':
    profile_dir = os.path.join(script_dir, 'Profiles')

### Check input
def is_fasta(input):
    with open(input, 'r') as handle:
        fa = SeqIO.parse(handle, 'fasta')
        return any(fa)

if (not is_fasta(fasta)) and check_inp:
    logging.critical('Input file is not in fasta format')
    sys.exit()

### Check output
out = os.path.join(out, '')
if not redo:
    try:
        os.mkdir(out)
    except FileExistsError:
        logging.error('Directory '+out+' already exists')
        sys.exit()

### Prodigal
if not aa:
    logging.info('Predicting ORFs with prodigal')
    with open(out + 'prodigal.out', 'w') as prodigal_out:
        with open(out + 'prodigal.log', 'w') as prodigal_log:
            subprocess.run(['prodigal', '-i', fasta, '-a', out+'proteins.faa', '-p', prod], stdout=prodigal_out, stderr=prodigal_log)

    # Check prodigal output
    if os.stat(out+'proteins.faa').st_size == 0:
        logging.critical('Prodigal failed. Check prodigal log')
        sys.exit()

    prot_path = out+'proteins.faa'

else:
    prot_path = fasta

### Clean up function
def clean(keep_prodigal):
    if not redo:
        logging.info('Removing temporary files')
        shutil.rmtree(out+'hmmer')
        os.remove(out+'hmmer.out')

        if os.stat(out+'hmmer.log').st_size == 0:
            os.remove(out+'hmmer.log')

        if not keep_prodigal:
            os.remove(out+'prodigal.out')
            os.remove(out+'proteins.faa')
            os.remove(out+'prodigal.log')

### Hmmer
if not redo:
    logging.info('Running HMMER against Cas profiles')
    os.mkdir(out+'hmmer')
    # Define function
    def hmmsearch(hmms, out):
        hmm_name = re.sub('\.hmm', '', hmms)
        with open(out+'hmmer.out', 'a') as hmmer_out:
            with open(out+'hmmer.log', 'a') as hmmer_log:
                subprocess.run(['hmmsearch', '--domtblout', os.path.join(out+'hmmer', hmm_name+'.tab'), os.path.join(profile_dir, hmms), prot_path], stdout=hmmer_out, stderr=hmmer_log)
    # Start multiprocess
    pool = mp.Pool(threads)
    # Each HMM
    hmm_dump = [pool.apply(hmmsearch, args=(hmms, out)) for hmms in os.listdir(profile_dir)]
    # Close multiprocess
    pool.close()

    # Combine
    logging.info('Parsing HMMER output')
    hmm_files = glob.glob(os.path.join(out+'hmmer', '*.tab'))
    with open(out+'hmmer.tab', 'w') as hmmer_tab:
        subprocess.run(['grep', '-v', '^#']+hmm_files, stdout=hmmer_tab)
        subprocess.run(['sed', '-i', 's/:/ /', out+'hmmer.tab'])

    hmm_df = pd.read_csv(out+'hmmer.tab', sep='\s+', header=None,
        usecols=(0,1,3,6,7,8,16,17,18,19,20,21,22,24,26),
        names=('Hmm','ORF','tlen','qlen','Eval','score','hmm_from','hmm_to','ali_from','ali_to','env_from','env_to','pprop','start','end'))

    hmm_df['Hmm'] = [re.sub('\.tab', '', re.sub(os.path.join(out, 'hmmer', ''), '', x)) for x in hmm_df['Hmm']]
    hmm_df.to_csv(out+'hmmer.tab', sep='\t', index=False)

else:
    hmm_df = pd.read_csv(out+'hmmer.tab', sep='\t')

if len(hmm_df) == 0:
    logging.info('No Cas proteins found. Exiting')
    clean(keep_prodigal)
    sys.exit()

# Add columns
# Acc
hmm_df['Acc'] = [re.sub("_[0-9]*$","",x) for x in hmm_df['ORF']]
# Gene position
hmm_df['Pos'] = [int(re.sub(".*_","",x)) for x in hmm_df['ORF']]
# Sequence coverage
hmm_df['Cov_seq'] = (hmm_df['ali_to'] - hmm_df['ali_from'] + 1) / hmm_df['tlen']
# Profile coverage
hmm_df['Cov_hmm'] = (hmm_df['hmm_to'] - hmm_df['hmm_from'] + 1) / hmm_df['qlen']

# Aggregate coverages of multiple aligments between similar HMMs and ORFs
hmm_df = hmm_df.groupby(['Hmm','ORF','tlen','qlen','Eval','score','start','end','Acc','Pos']).agg({'Cov_seq':'sum','Cov_hmm':'sum'}).reset_index()

# Pick best hit
hmm_df.sort_values('score', ascending=False, inplace=True)
hmm_df.drop_duplicates('ORF', inplace=True)

# Typing
logging.info('Subtyping putative operons')
out_df = typing(dat=hmm_df,
    out=out,
    threads=threads,
    scoring=scoring,
    dist=dist,
    overall_eval=oev,
    overall_cov_seq=ocs,
    overall_cov_hmm=och,
    one_gene_eval=sev,
    one_gene_cov_seq=scs,
    one_gene_cov_hmm=sch,
    two_gene_eval=tev,
    two_gene_cov_seq=tcs,
    two_gene_cov_hmm=tch,
    single_gene_types=('II-A','II-B','II-C','V-A','V-B','V-C','V-D','V-E','V-F','V-G','V-H','V-I','V-J','VI-A','VI-B1','VI-B2','VI-C','VI-D'),
    VF_eval=vfe,
    VF_cov_hmm=vfc)

if len(out_df) == 0:
    logging.info('No operons found. Exiting')
    clean(keep_prodigal)
    sys.exit()

# Remove temporary hmmer files
clean(keep_prodigal)

# Output
logging.info('Done')
operons_good = out_df[~out_df['Prediction'].isin(['False', 'Ambiguous', 'Partial'])]
operons_put = out_df[out_df['Prediction'].isin(['False', 'Ambiguous', 'Partial'])]

operons_good.to_csv(out+'cas_operons.tab', sep='\t', index=False)
operons_put.to_csv(out+'cas_operons_putative.tab', sep='\t', index=False)
