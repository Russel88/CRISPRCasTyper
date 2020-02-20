import os
import logging
import sys
import shutil

import pandas as pd

from Bio import SeqIO

class Controller(object):

    def __init__(self, args):
       
        self.fasta = args.input
        self.out = args.output
        self.threads = args.threads
        self.dist = args.dist
        self.prod = args.prodigal
        self.scoring = args.scores
        self.pdir=args.hmms
        self.oev = args.overall_eval
        self.ocs = args.overall_cov_seq
        self.och = args.overall_cov_hmm
        self.tev = args.two_gene_eval
        self.tcs = args.two_gene_cov_seq
        self.tch = args.two_gene_cov_hmm
        self.sev = args.single_gene_eval
        self.scs = args.single_gene_cov_seq
        self.sch = args.single_cov_hmm
        self.vfe = args.vf_eval
        self.vfc = args.vf_cov_hmm
        self.check_inp = args.skip_check
        self.keep_prodigal = args.keep_prodigal
        self.lvl = args.log_lvl
        self.aa = args.aa
        self.redo = args.redo_typing

        # Logger
        logging.basicConfig(format='\033[36m'+'[%(asctime)s] %(levelname)s:'+'\033[0m'+' %(message)s', datefmt='%Y-%m-%d %H:%M:%S', level=self.lvl)
        logging.info('Running CasPredict version 0.2.1')

        # Force consistency
        self.out = os.path.join(self.out, '')

        if self.redo:
            self.aa = True
            self.check_inp = True

        if self.aa:
            self.keep_prodigal = True
            self.prot_path = self.fasta
        else:
            self.prot_path = self.out+'proteins.faa'

        # Hardcode single gene types
        self.single_gene_types=('II-A','II-B','II-C','V-A','V-B','V-C','V-D','V-E','V-F','V-G','V-H','V-I','V-J','VI-A','VI-B','VI-C','VI-D')

        # Check databases
        self.check_db()
        
        # Check input and output
        self.check_input()
        self.check_out()

    def check_out(self):

        if not self.redo:
            try:
                os.mkdir(self.out)
            except FileExistsError:
                logging.error('Directory '+self.out+' already exists')
                sys.exit()

    def check_input(self):

        if not self.check_inp:
            if os.path.isfile(self.fasta):
                if not self.is_fasta():
                    logging.error('Input file is not in fasta format')
                    sys.exit()
            else:
                logging.error('Could not find input file')
                sys.exit()

    def is_fasta(self):
        with open(self.fasta, 'r') as handle:
            fa = SeqIO.parse(handle, 'fasta')
            return any(fa)

    def clean(self):
        if not self.redo:
            logging.info('Removing temporary files')
            shutil.rmtree(self.out+'hmmer')

            if os.stat(self.out+'hmmer.log').st_size == 0:
                os.remove(self.out+'hmmer.log')

            if not self.keep_prodigal:
                os.remove(self.out+'proteins.faa')
                os.remove(self.out+'prodigal.log')
    
    def check_db(self):
       # Try to load CasScoring table
        if os.path.isfile(self.scoring):
            try:
                dump = pd.read_csv(self.scoring, sep=",")
            except:
                logging.error('CasScoring table could not be loaded')
                sys.exit()
        else:
            logging.error('CasScoring table could not be found')
            sys.exit()

        # Look if HMM profile dir exists
        if os.path.isdir(self.pdir):
            for i in os.listdir(self.pdir):
                if not i.lower().endswith('.hmm'):
                    logging.error('There are non-HMM profiles in the HMM profile directory')
                    sys.exit()
        else:
            logging.error('Could not find HMM profile directory')
            sys.exit()
