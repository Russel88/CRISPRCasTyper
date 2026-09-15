from __future__ import annotations

import os
import glob
import logging
import sys
import shutil
import json
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
import pandas as pd

from Bio import SeqIO
from cctyper.resources import resolve_database_path

class Controller(object):

    def __init__(self, args):
       
        self.fasta: str = args.input
        self.out: str = args.output
        self.threads = args.threads
        self.dist = args.dist
        self.gff = args.gff
        self.prot = args.prot
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
        self.repeat_search_backend = args.repeat_search_backend

        # MinCED
        self.searchWL = args.searchWL
        self.minNR = args.minNR
        self.minRL = args.minRL
        self.maxRL = args.maxRL
        self.minSL = args.minSL
        self.maxSL = args.maxSL

        self.any_cas = False
        self.any_operon = False
        self.any_crispr = False

        # Logger
        try:
            cctyper_version = version("cctyper")
        except PackageNotFoundError:
            cctyper_version = "unknown"

        if self.simplelog:
            logging.basicConfig(format='[%(asctime)s] %(levelname)s: %(message)s', datefmt='%Y-%m-%d %H:%M:%S', level=self.lvl)
        else:
            logging.basicConfig(format='\033[36m'+'[%(asctime)s] %(levelname)s:'+'\033[0m'+' %(message)s', datefmt='%Y-%m-%d %H:%M:%S', level=self.lvl)
        logging.info('Running CRISPRCasTyper version {}'.format(cctyper_version))

        # Pick repeat search backend now that logging is up
        self.resolve_backend(args)

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

    def resolve_backend(self, args):
        '''
        Prefer blast when installed, else fall back to the bundled edit distance backend
        '''

        has_blast = shutil.which('blastn') is not None and shutil.which('makeblastdb') is not None

        if self.repeat_search_backend == 'auto':
            self.repeat_search_backend = 'blast' if has_blast else 'myers'
            logging.debug('Repeat search backend: %s', self.repeat_search_backend)
        elif self.repeat_search_backend == 'blast' and not has_blast:
            logging.error('--repeat_search_backend blast requires blastn and makeblastdb in PATH')
            sys.exit(1)

        # Record what actually ran, not 'auto', so arguments.tab is reproducible
        args.repeat_search_backend = self.repeat_search_backend

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

        if self.gff or self.prot:
            if not (self.gff and self.prot):
                logging.error('Both --gff and --prot must be provided together')
                sys.exit()
            if not os.path.isfile(self.gff):
                logging.error('Could not find GFF file %s', self.gff)
                sys.exit()
            if not os.path.isfile(self.prot):
                logging.error('Could not find protein FASTA file %s', self.prot)
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
            fixed_path = Path(self.out) / 'fixed_input.fna'
            with open(self.fasta, 'r') as src, fixed_path.open('w') as dst:
                for line in src:
                    if line.startswith('>'):
                        dst.write('>Contig'+line[1:])
                    else:
                        dst.write(line)
            self.fasta = str(fixed_path)
            self.len_dict = {'Contig'+str(key): val for key, val in self.len_dict.items()}
            self.seq_dict = {'Contig'+str(key): val for key, val in self.seq_dict.items()}

    def clean(self):
        if not self.redo:

            if self.num_headers:
                os.remove(self.out+'fixed_input.fna')

            if os.path.exists(self.out+'hmmer.log') and os.stat(self.out+'hmmer.log').st_size == 0:
                os.remove(self.out+'hmmer.log')

            if self.customhmm != '':
                if os.stat(self.out+'hmmer_custom.log').st_size == 0:
                    os.remove(self.out+'hmmer_custom.log')
                
            if not self.keep_tmp:
                
                logging.info('Removing temporary files')
                
                shutil.rmtree(self.out+'hmmer', ignore_errors=True)

                for tmp_file in ['minced.out', 'prodigal.log', 'proteins.faa']:
                    tmp_path = self.out + tmp_file
                    if os.path.exists(tmp_path):
                        os.remove(tmp_path)

                if os.path.exists(self.out+'blast.tab'):
                    os.remove(self.out+'blast.tab')
                    for fl in glob.glob(self.out+'Flank.*'):
                        os.remove(fl)

    def check_db(self):

        try:
            self.db = resolve_database_path(self.db, logger=logging.getLogger(__name__))
        except RuntimeError as err:
            logging.error(err)
            sys.exit()

        self.scoring = os.path.join(self.db, 'CasScoring.csv')
        self.hmm_compressed = os.path.join(self.db, 'cctyper_profiles.hmm.zst')
        self.hmm_plain = os.path.join(self.db, 'cctyper_profiles.hmm')
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

        # Locate HMM profile database (zstd preferred)
        if os.path.isfile(self.hmm_plain):
            self.hmm_db = self.hmm_plain
        elif os.path.isfile(self.hmm_compressed):
            self.hmm_db = self.hmm_compressed
        else:
            logging.error('Could not find HMM profile database (expected %s or %s)', self.hmm_compressed, self.hmm_plain)
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
