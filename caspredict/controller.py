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
        self.db = args.db
        self.circular = args.circular
        self.oev = args.overall_eval
        self.ocs = args.overall_cov_seq
        self.och = args.overall_cov_hmm
        self.check_inp = args.skip_check
        self.keep_tmp = args.keep_tmp
        self.lvl = args.log_lvl
        self.redo = args.redo_typing
        self.kmer = args.kmer
        self.crispr_cas_dist = args.ccd
        self.pred_prob = args.pred_prob
        self.noplot = args.no_plot
        self.scale = args.scale
        self.nogrid = args.no_grid
        self.expand = args.expand
        self.plotexpand = args.plot_expand
        self.simplelog = args.simplelog

        self.any_cas = False
        self.any_operon = False
        self.any_crispr = False

        # Logger
        if self.simplelog:
            logging.basicConfig(format='[%(asctime)s] %(levelname)s: %(message)s', datefmt='%Y-%m-%d %H:%M:%S', level=self.lvl)
        else:
            logging.basicConfig(format='\033[36m'+'[%(asctime)s] %(levelname)s:'+'\033[0m'+' %(message)s', datefmt='%Y-%m-%d %H:%M:%S', level=self.lvl)
        logging.info('Running CasPredict version 0.5.4')

        # Force consistency
        self.out = os.path.join(self.out, '')

        if self.redo:
            self.check_inp = True

        self.prot_path = self.out+'proteins.faa'

        # Check databases
        self.check_db()
        
        # Check input and output
        self.check_input()
        self.check_out()

        # If redo check if any crisprs and operons
        if self.redo:
            if os.path.exists(self.out+'cas_operons.tab') or os.path.exists(self.out+'cas_operons_putative.tab'):
                self.any_operon = True
            if os.path.exists(self.out+'crisprs_all.tab'):
                self.any_crispr = True

        # Write arguments
        da = vars(args)
        f = open(self.out+'arguments.tab', 'w')
        for k, v in da.items():
            f.write('{}:\t{}\n'.format(k, v))
        f.close()

        # If circular get lengths
        if self.circular:
            self.get_length()

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
        
        try:
            with open(self.fasta, 'r') as handle:
                fa = SeqIO.parse(handle, 'fasta')
                [float(x.id) for x in fa]
                logging.error('Numeric fasta headers not supported')
                return False
        except:
            with open(self.fasta, 'r') as handle:
                fa = SeqIO.parse(handle, 'fasta')
                return any(fa)

    def clean(self):
        if not self.redo:

            if os.stat(self.out+'hmmer.log').st_size == 0:
                os.remove(self.out+'hmmer.log')

            if not self.keep_tmp:
                
                logging.info('Removing temporary files')
                
                shutil.rmtree(self.out+'hmmer')
                
                os.remove(self.out+'minced.out')
                os.remove(self.out+'prodigal.log')
                os.remove(self.out+'proteins.faa')

    def check_db(self):
        
        if self.db == '':
            try:
                self.db = os.environ['CASPREDICT_DB']
            except:
                logging.error('Could not find database directory')
                sys.exit()

        self.scoring = os.path.join(self.db, 'CasScoring.csv')
        self.pdir = os.path.join(self.db, 'Profiles', '')
        self.xgb = os.path.join(self.db, "xgb_repeats.model")
        self.typedict = os.path.join(self.db, "type_dict.tab")
        self.cutoffdb = os.path.join(self.db, "cutoffs.tab")

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

        # Load specific cutoffs
        with open(self.cutoffdb, 'r') as f:
            rs = (ll.rstrip().split(':') for ll in f)
            self.cutoffs = {r[0].lower():r[1].split(',') for r in rs}

    def get_length(self):
        with open(self.fasta, 'r') as handle:
            self.len_dict = {}
            for fa in SeqIO.parse(handle, 'fasta'):
                self.len_dict[fa.id] = len(fa.seq)
        
