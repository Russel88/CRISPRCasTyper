import os
import subprocess
import logging
import sys
import re
import math
import random

import statistics as st

from Bio import pairwise2
from joblib import Parallel, delayed

# Define the CRISPR class
class CRISPR(object):
    count = 0
    def __init__(self, sequence, exact_stats):
        self.sequence = sequence.rstrip()
        CRISPR.count += 1
        self.crispr = '{}_{}'.format(self.sequence, CRISPR.count)
        self.repeats = []
        self.spacers = []
        self.exact = exact_stats
    def setPos(self, start, end):
        self.start = start.rstrip()
        self.end = end.rstrip()
    def addRepeat(self, repeat):
        self.repeats.append(repeat.rstrip())
    def addSpacer(self, spacer):
        self.spacers.append(spacer.rstrip())
    def getConsensus(self):
        self.cons = max(set(self.repeats), key = self.repeats.count) 
    def identity(self, i, j, sqlst):
        align = pairwise2.align.globalxx(sqlst[i], sqlst[j])
        return(align[0][2]/align[0][4]*100)
    def identLoop(self, seqs, threads):
        if self.exact:
            sqr = range(len(seqs))
        else:
            if len(seqs) > 10:
                n_samp = 10
            else:
                n_samp = len(seqs)
            sqr = random.sample(range(len(seqs)), n_samp)
        idents = Parallel(n_jobs=threads)(delayed(self.identity)(k, l, seqs) for k in sqr for l in sqr if k > l)
        return(st.mean(idents))
    def stats(self, threads, rep_id, spa_id, spa_sem):
        self.spacer_identity = round(self.identLoop(self.spacers, threads), 1)
        self.repeat_identity = round(self.identLoop(self.repeats, threads), 1)
        self.spacer_len = round(st.mean([len(x) for x in self.spacers]), 1)
        self.repeat_len = round(st.mean([len(x) for x in self.repeats]), 1)
        self.spacer_sem = round(st.stdev([len(x) for x in self.spacers])/math.sqrt(len(self.spacers)), 1)
        self.trusted = (self.repeat_identity > rep_id) & (self.spacer_identity < spa_id) & (self.spacer_sem < spa_sem)

class Minced(object):
    
    def __init__(self, obj):
        self.master = obj
        for key, val in vars(obj).items():
            setattr(self, key, val)

    def run_minced(self):

        if not self.redo:
            logging.info('Predicting CRISPR arrays with minced')

            # Run minced
            subprocess.run(['minced', 
                            self.fasta, 
                            self.out+'minced.out'], 
                            stdout=subprocess.DEVNULL, 
                            stderr=subprocess.DEVNULL)
        
            # Parse
            self.parse_minced()

            # Write results
            self.write_crisprs()
            self.write_spacers()

    def parse_minced(self):
        file = open(self.out+'minced.out', 'r')

        random.seed(self.seed)
        
        crisprs = []
        for ll in file:
            # Record sequence accession
            if ll.startswith('Sequence'):
                sequence_current = re.sub('\' \(.*', '', re.sub('Sequence \'', '', ll))
            # Create instance of CRISPR and add positions
            if ll.startswith('CRISPR'):
                crisp_tmp = CRISPR(sequence_current, self.exact_stats)
                pos = re.sub('.*Range: ', '', ll)
                start = re.sub(' - .*', '', pos)
                end = re.sub('.* - ', '', pos)
                crisp_tmp.setPos(start, end)
            # Add Repeats and Spacers to the current instance
            if ll[:1].isdigit():
                lll = ll.split()
                if len(lll) == 7:
                    crisp_tmp.addRepeat(lll[1])
                    crisp_tmp.addSpacer(lll[2])
                if len(lll) == 2:
                    crisp_tmp.addRepeat(lll[1])
            # Save the instance
            if ll.startswith('Repeats'):
                crisp_tmp.getConsensus()
                crisp_tmp.stats(self.threads, self.repeat_id, self.spacer_id, self.spacer_sem)
                crisprs.append(crisp_tmp)

        file.close()

        self.crisprs = crisprs

    def write_crisprs(self):
       
        header = '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format('Contig',
                                                                    'CRISPR',
                                                                    'Start',
                                                                    'End',
                                                                    'Consensus_repeat',
                                                                    'N_repeats',
                                                                    'Repeat_len',
                                                                    'Spacer_len_avg',
                                                                    'Repeat_identity',
                                                                    'Spacer_identity',
                                                                    'Spacer_len_sem',
                                                                    'Trusted')
       
        def write_crisp(handle, cris):
            handle.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(cris.sequence,
                                                   cris.crispr,
                                                   cris.start,
                                                   cris.end,
                                                   cris.cons,
                                                   len(cris.repeats),
                                                   cris.repeat_len,
                                                   cris.spacer_len,
                                                   cris.repeat_identity,
                                                   cris.spacer_identity,
                                                   cris.spacer_sem,
                                                   cris.trusted))
            
        f = open(self.out+'crisprs_all.tab', 'w')
        f.write(header)
        for crisp in self.crisprs:
            write_crisp(f, crisp)
        f.close()

    def write_spacers(self):
        
        if len(self.crisprs) > 0:
            os.mkdir(self.out+'spacers')
            for crisp in self.crisprs:
                f = open(self.out+'spacers/{}.fa'.format(crisp.crispr), 'w')
                n = 0
                for sq in crisp.spacers:
                    n += 1
                    f.write('>{}:{}\n'.format(crisp.crispr, n))
                    f.write('{}\n'.format(sq))

                f.close()




