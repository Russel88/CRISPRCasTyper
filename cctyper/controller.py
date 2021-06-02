import os
import logging
import sys
import shutil
import json
import pkg_resources
import subprocess

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
        self.keep_tmp = args.keep_tmp
        self.lvl = args.log_lvl
        self.redo = args.redo_typing
        self.kmer = args.kmer
        self.crispr_cas_dist = args.ccd
        self.pred_prob = args.pred_prob
        self.noplot = args.no_plot
        self.nogrid = args.no_grid
        self.expand = args.expand
        self.simplelog = args.simplelog
        self.customhmm = args.custom_hmm
        self.repeat_id = args.repeat_id
        self.spacer_id = args.spacer_id
        self.spacer_sem = args.spacer_sem
        self.exact_stats = args.exact_stats
        self.seed = args.seed
        self.skip_blast = args.skip_blast

        self.any_cas = False
        self.any_operon = False
        self.any_crispr = False

        # Logger
        if self.simplelog:
            logging.basicConfig(format='[%(asctime)s] %(levelname)s: %(message)s', datefmt='%Y-%m-%d %H:%M:%S', level=self.lvl)
        else:
            logging.basicConfig(format='\033[36m'+'[%(asctime)s] %(levelname)s:'+'\033[0m'+' %(message)s', datefmt='%Y-%m-%d %H:%M:%S', level=self.lvl)
        logging.info('Running CRISPRCasTyper version {}'.format(pkg_resources.require("cctyper")[0].version))

        # kmer warning
        if self.kmer != 4:
            logging.warning('kmer argument should only be used if the repeatTyper model is trained with a different kmer than 4.')
        
        # Force consistency
        self.out = os.path.join(self.out, '')

        self.prot_path = self.out+'proteins.faa'

        # Check databases
        self.check_db()
        
        # Check input and output
        self.check_out()
        self.check_input()

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

    def check_out(self):

        if not self.redo:
            try:
                os.mkdir(self.out)
            except FileExistsError:
                logging.error('Directory '+self.out+' already exists')
                sys.exit()

    def check_input(self):

        if os.path.isfile(self.fasta):
            self.check_fasta()
        else:
            logging.error('Could not find input file')
            sys.exit()

    def check_fasta(self):
        
        # Get sequence lengths
        with open(self.fasta, 'r') as handle:
            self.len_dict = {}
            self.seq_dict = {}
            for fa in SeqIO.parse(handle, 'fasta'):
                if fa.id in self.len_dict:
                    logging.error('Duplicate fasta headers detected')
                    sys.exit()
                self.len_dict[fa.id] = len(fa.seq)
                self.seq_dict[fa.id] = fa.seq
            
        # Check for numeric headers
        self.num_headers = False
        for i in self.len_dict.keys():
            try:
                dump = float(i)
                self.num_headers = True
            except:
                pass
        
        if self.num_headers:
            logging.warning('Numeric fasta headers detected. A prefix is added to the names')
            new_fasta = open(self.out+'fixed_input.fna', 'w')
            subprocess.run(['sed', 's/^>/>Contig/', self.fasta], stdout = new_fasta)
            new_fasta.close()
            self.fasta = self.out+'fixed_input.fna'
            self.len_dict = {'Contig'+str(key): val for key, val in self.len_dict.items()}
            self.seq_dict = {'Contig'+str(key): val for key, val in self.seq_dict.items()}

    def clean(self):
        if not self.redo:

            if self.num_headers:
                os.remove(self.out+'fixed_input.fna')

            if os.stat(self.out+'hmmer.log').st_size == 0:
                os.remove(self.out+'hmmer.log')

            if self.customhmm != '':
                if os.stat(self.out+'hmmer_custom.log').st_size == 0:
                    os.remove(self.out+'hmmer_custom.log')
                
            if not self.keep_tmp:
                
                logging.info('Removing temporary files')
                
                shutil.rmtree(self.out+'hmmer')
                
                os.remove(self.out+'minced.out')
                os.remove(self.out+'prodigal.log')
                os.remove(self.out+'proteins.faa')

                if os.path.exists(self.out+'blast.tab'):
                    os.remove(self.out+'blast.tab')
                    os.remove(self.out+'Flank.fna')
                    os.remove(self.out+'Flank.nhr')
                    os.remove(self.out+'Flank.nin')
                    os.remove(self.out+'Flank.nsq')

    def check_db(self):
        
        if self.db == '':
            try:
                self.db = os.environ['CCTYPER_DB']
            except:
                logging.error('Could not find database directory')
                sys.exit()

        self.scoring = os.path.join(self.db, 'CasScoring.csv')
        self.pdir = os.path.join(self.db, 'Profiles', '')
        self.xgb = os.path.join(self.db, "xgb_repeats.model")
        self.typedict = os.path.join(self.db, "type_dict.tab")
        self.cutoffdb = os.path.join(self.db, "cutoffs.tab")
        self.ifdb = os.path.join(self.db, "interference.json")
        self.addb = os.path.join(self.db, "adaptation.json")
        self.repeatdb = os.path.join(self.db, "repeats.fa")

        # Load CasScoring table
        if os.path.isfile(self.scoring):
            try:
                self.scores = pd.read_csv(self.scoring, sep=",")
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

        # Load mandatory gene files
        with open(self.ifdb, 'r') as f:
            self.compl_interf = json.load(f)
        with open(self.addb, 'r') as f:
            self.compl_adapt = json.load(f)

