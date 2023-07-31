import os
import subprocess
import sys
import logging
import re
import glob
import tqdm

import multiprocess as mp
import pandas as pd

class HMMER(object):
    
    def __init__(self, obj):
        self.master = obj
        for key, val in vars(obj).items():
            setattr(self, key, val)

    def main_hmm(self):
        # If redo just get the table
        if self.redo:
            self.read_hmm()
        # Else run HMMER load and write data
        else:
            self.run_hmm()
            self.load_hmm()
            self.write_hmm()
            self.run_custom_hmm()

        # Check if any cas genes
        self.check_hmm()

        # Parse
        self.parse_hmm()

        # Load Custom HMM db
        self.load_custom_hmm()

    # A single search
    def hmmsearch(self, hmm):
       
        hmm_name = re.sub('\.hmm', '', hmm)
        
        logging.debug('Running HMMER against '+hmm_name)
        
        with open(self.out+'hmmer.log', 'a') as hmmer_log:
            subprocess.run(['hmmsearch', 
                            '--domtblout', os.path.join(self.out+'hmmer', hmm_name+'.tab'), 
                            os.path.join(self.pdir, hmm), 
                            self.prot_path], 
                            stdout=subprocess.DEVNULL, 
                            stderr=hmmer_log)

    # Parallel search of all HMMs
    def run_hmm(self):
        
        logging.info('Running HMMER against Cas profiles')
        
        # Make dir
        os.mkdir(self.out+'hmmer')
        # Start multiprocess
        pool = mp.Pool(self.threads)
        # Each HMM
        if self.lvl == 'DEBUG' or self.simplelog:
            list(pool.imap(self.hmmsearch, os.listdir(self.pdir)))
        else:
            list(tqdm.tqdm(pool.imap(self.hmmsearch, os.listdir(self.pdir)), total=len(os.listdir(self.pdir))))
        # Close multiprocess
        pool.close()

    # Load data
    def load_hmm(self):
    
        logging.debug('Loading HMMER output')
        
        # Get files
        hmm_files = glob.glob(os.path.join(self.out+'hmmer', '*.tab'))
        # Parse externally
        with open(self.out+'hmmer.tab', 'w') as hmmer_tab:
            subprocess.run(['grep', '-v', '^#']+hmm_files, stdout=hmmer_tab)
            subprocess.run(['sed', '-i', 's/:/ /', self.out+'hmmer.tab'])

        # Load
        usecols = [0,1,3,6,7,
                    8,16,17,18,19,
                    20,21,22]
        names = ['Hmm','ORF','tlen','qlen','Eval',
                'score','hmm_from','hmm_to','ali_from','ali_to',
                'env_from','env_to','pprop']
        if not (self.gff and self.prot):
            usecols = usecols + [24,26,28]
            names= names + ['start','end','strand'] 
        print(names)
        hmm_df = pd.read_csv(self.out+'hmmer.tab', sep='\s+', header=None,
            usecols=usecols,
            names=names)
        print(hmm_df)
        # Parse HMM names
        hmm_df['Hmm'] = [re.sub('\.tab', '', 
                        re.sub(os.path.join(self.out, 'hmmer', ''), '', x)) 
                        for x in hmm_df['Hmm']]
        
        if self.gff and self.prot:
            self.genes = pd.read_csv(self.out+'genes.tab', sep='\t')
            hmm_df = pd.merge(hmm_df,self.genes[["Start","End","Strand","Contig","Pos","protein_id"]],
                                            left_on="ORF",right_on="protein_id",how="left").drop("protein_id",axis=1)
            hmm_df.rename(columns={'Contig': 'Acc','Strand': 'strand', 'Start': 'start', 'End': 'end'}, inplace=True)
        else:
            hmm_df['Acc'] = [re.sub("_[0-9]*$","",x) for x in hmm_df['ORF']]
            hmm_df['Pos'] = [int(re.sub(".*_","",x)) for x in hmm_df['ORF']]

        # Coverages of aligments
        def covs(df_sub):
            df_sub['Cov_seq'] = len(set([x for sublst in [list(range(i,j)) 
                for i,j in zip(df_sub['ali_from'], df_sub['ali_to']+1)] 
                for x in sublst])) / df_sub['tlen']
            df_sub['Cov_hmm'] = len(set([x for sublst in [list(range(i,j)) 
                for i,j in zip(df_sub['hmm_from'], df_sub['hmm_to']+1)] 
                for x in sublst])) / df_sub['qlen']
            df_sub = df_sub[['Hmm','ORF','tlen','qlen','Eval','score',
                            'start','end','Acc','Pos','Cov_seq','Cov_hmm','strand']]
            df_sub = df_sub.drop_duplicates()
            return df_sub

        hmm_df = hmm_df.groupby(['Hmm','ORF']).apply(covs)
        hmm_df.reset_index(drop=True, inplace=True)
        self.hmm_df = hmm_df.drop_duplicates()
    # Write to file
    def write_hmm(self):
        self.hmm_df.to_csv(self.out+'hmmer.tab', sep='\t', index=False)

    # Read from file
    def read_hmm(self):
        try:        
            self.hmm_df = pd.read_csv(self.out+'hmmer.tab', sep='\t')
        except:
            logging.error('No matches to Cas HMMs')
            sys.exit()

    # Check if any cas genes
    def check_hmm(self):
        if len(self.hmm_df) == 0:
            logging.info('No Cas proteins found.')
        else:
            self.any_cas = True

    # Parse
    def parse_hmm(self):

        if self.any_cas:
        
            logging.debug('Parsing HMMER output')

            # Pick best hit
            self.hmm_df.sort_values('score', ascending=False, inplace=True)
            self.hmm_df.drop_duplicates('ORF', inplace=True)

    def run_custom_hmm(self):
        
        if self.customhmm != '':
            logging.info('Running HMMER against custom HMM profiles')
            
            with open(self.out+'hmmer_custom.log', 'a') as hmmer_log:
                subprocess.run(['hmmsearch', 
                                '--tblout', self.out+'hmmer_custom.tab', 
                                '--cpu', str(self.threads),
                                self.customhmm, 
                                self.prot_path], 
                                stdout=subprocess.DEVNULL, 
                                stderr=hmmer_log)
            
    def load_custom_hmm(self):

        if self.customhmm != '':

            # Check if successful
            if not os.path.isfile(self.out+'hmmer_custom.tab'):
                logging.error('HMMER failed running on the custom HMM database')
                sys.exit()

            # Load
            self.custom_hmm_df = pd.read_csv(self.out+'hmmer_custom.tab', sep='\s+', comment='#', 
                header=None, usecols=(0, 2, 3, 4, 5), 
                names=('Target', 'Query', 'Acc', 'E-value', 'Score'))
                        
            # Remove low E-value hits
            self.custom_hmm_df = self.custom_hmm_df[self.custom_hmm_df['E-value'] < self.oev]
            
            # Pick best hit
            self.custom_hmm_df.sort_values('Score', ascending=False, inplace=True)
            self.custom_hmm_df.drop_duplicates('Target', inplace=True)

            # New columns
            if self.gff and self.prot:
                self.genes = pd.read_csv(self.out+'genes.tab', sep='\t')
                self.custom_hmm_df = pd.merge(self.custom_hmm_df,self.genes[["Contig","Pos","protein_id"]],
                                              left_on="Target",right_on="protein_id",how="left").drop("protein_id",axis=1)
            else:
                self.custom_hmm_df['Contig'] = [re.sub("_[0-9]*$","",x) for x in self.custom_hmm_df['Target']]
                self.custom_hmm_df['Pos'] = [int(re.sub(".*_","",x)) for x in self.custom_hmm_df['Target']]

