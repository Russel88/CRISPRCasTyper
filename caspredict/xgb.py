import logging
import os

import pandas as pd
import numpy as np
import itertools as it
import xgboost as xgb

class XGB(object):

    def __init__(self, obj):
        self.master = obj
        for key, val in vars(obj).items():
            setattr(self, key, val)
        
        base_for = "ACGT"
        base_rev = "TGCA"
        self.comp_tab = str.maketrans(base_for, base_rev)
    
    def load_xgb_model(self):

        logging.debug('Loading xgboost model')

        bst = xgb.Booster({'nthread':self.threads})
        bst.load_model(self.xgb)
        self.bst = bst

        # Load label dict here:
        with open(self.typedict, 'r') as f:
            rs = (ll.rstrip().split(':') for ll in f)
            self.label_dict = {r[1]:r[0] for r in rs} 

    def generate_canonical_kmer(self):

        logging.debug('Generating canonical {}mers'.format(self.kmer))

        letters = ['A','C','G','T']
        all_kmer = [''.join(k) for k in it.product(letters, repeat=self.kmer)]
        all_kmer_rev = [x.translate(self.comp_tab)[::-1] for x in all_kmer]
        can_kmer = list(it.compress(all_kmer_rev, [not kf < kr for kf,kr in zip(all_kmer, all_kmer_rev)]))
        can_kmer.sort()
        self.can_kmer = can_kmer

    def count_kmer(self, seq):
        kmer_d = {}
        for i in range(len(seq) - self.kmer + 1):
            kmer_for = seq[i:(i+self.kmer)]
            kmer_rev = kmer_for.translate(self.comp_tab)[::-1]
            if kmer_for < kmer_rev: 
                kmer = kmer_for
            else: 
                kmer = kmer_rev
            if kmer in kmer_d:
                kmer_d[kmer] += 1
            else: 
                kmer_d[kmer] = 1
        return kmer_d

    def xgb_run(self):
        
        if not self.redo:
            logging.info('Predicting subtype of CRISPR repeats')
            
            # Get repeats
            self.repeats = [x.cons for x in self.crisprs]

            # Predict
            self.predict_repeats()

            # Add to file
            df = pd.read_csv(self.out+'crisprs_all.tab', sep='\t')
            df['Prediction'] = self.z_type
            df['Subtype'] = self.z_type
            df['Subtype_probability'] = self.z_max
            df.loc[df.Subtype_probability < self.pred_prob, 'Prediction'] = 'Unknown'
            df['Subtype_probability'] = df['Subtype_probability'].round(3)
           
            if len(df) > 0:
                df.to_csv(self.out+'crisprs_all.tab', sep='\t', index=False)
                self.any_crispr = True
            else:
                logging.info('No CRISPRs found.')
                os.remove(self.out+'crisprs_all.tab')
    
    def predict_repeats(self):

        # Prepare
        self.load_xgb_model()
        self.generate_canonical_kmer()
        self.repeats = [x.upper() for x in self.repeats]

        # Count kmers (first index is a to ensure all kmers are in the df)
        z_df = pd.DataFrame([dict(zip(self.can_kmer, np.zeros(len(self.can_kmer))))] + [self.count_kmer(x) for x in self.repeats]).fillna(0)
        z_df = z_df.reindex(sorted(z_df.columns), axis=1)
        
        # Predict
        self.z_pred = self.bst.predict(xgb.DMatrix(z_df), ntree_limit=int(self.bst.attr('best_iteration')))
        
        # Get type and max probability
        self.z_best = [x.argmax() for x in self.z_pred][1:len(self.z_pred)]
        self.z_max = [x.max() for x in self.z_pred][1:len(self.z_pred)]

        # Convert to type string
        self.z_type = [self.label_dict[str(x)] for x in self.z_best]

    def print_xgb(self):
        
        for i in range(len(self.repeats)):
            print('{}\t{}\t{}'.format(self.repeats[i], 
                                    self.z_type[i],
                                    self.z_max[i]))


