import os
import sys
import re

import pandas as pd
import itertools as it
import numpy as np
import xgboost as xgb

from sklearn.model_selection import train_test_split
from sklearn.metrics import confusion_matrix

from cctyper.xgb import XGB

class XGBTrain(object):

    def __init__(self, args):
       
        self.input = args.input
        self.out = args.output
        self.kmer = args.kmer
        self.minr = args.minr
        self.rnd_seed = args.rnd_seed
        self.test_size = args.test_size
        self.eta = args.eta
        self.threads = args.threads
        self.num_rounds = args.num_rounds
        self.early_stop = args.early_stop
        self.max_depth = args.max_depth
        self.subsample = args.subsample
        self.colsample_bytree = args.colsample_bytree
        self.nfold = args.nfold

        # Alphabet
        base_for = "ACGT"
        base_rev = "TGCA"
        self.comp_tab = str.maketrans(base_for, base_rev)
        
        # Read input
        self.read_input()

        # Output
        self.out = os.path.join(self.out, '')
        self.check_out()

        # Prune tybes
        self.prune_input()

        # Prepare
        self.prepare_data()

        # Train
        self.train()

        # Test
        self.test()

    def read_input(self):
        
        # Load input:i
        self.dat = pd.read_csv(self.input, header=None, sep='\t', names=('Type',"Seq"))

        # Check input
        def is_dna(s):
            match = re.match("^[ACTGactg]*$", s)
            return match is not None

        for rep in list(self.dat['Seq']):
            if not is_dna(rep):
                print('Error - Non-DNA letters found in sequence:')
                print(rep)
                sys.exit()

    def check_out(self):

        try:
            os.mkdir(self.out)
        except FileExistsError:
            print('Directory '+self.out+' already exists')
            sys.exit()

    def prune_input(self):
        
        subtype, count = np.unique(self.dat['Type'], return_counts=True)

        print('\033[92m' + 'Counts of subtypes:' + '\033[0m')
        print(list(zip(subtype, count)))

        self.incl = list(subtype[count >= self.minr])

        print('\033[92m' + "Including the following subtypes:" + '\033[0m')
        print(self.incl)

        # Dictionary
        self.label_dict = dict(zip(self.incl, range(len(self.incl))))
        
        f = open(self.out+'type_dict.tab', 'w')
        for k, v in self.label_dict.items():
            f.write('{}:{}\n'.format(k, v))
        f.close()

        self.dat = self.dat[self.dat['Type'].isin(self.incl)]

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

    def generate_canonical_kmer(self):

        letters = ['A','C','G','T']
        all_kmer = [''.join(k) for k in it.product(letters, repeat=self.kmer)]
        all_kmer_rev = [x.translate(self.comp_tab)[::-1] for x in all_kmer]
        can_kmer = list(it.compress(all_kmer_rev, [not kf < kr for kf,kr in zip(all_kmer, all_kmer_rev)]))
        can_kmer.sort()
        self.can_kmer = can_kmer
    
    def prepare_data(self):

        self.generate_canonical_kmer()

        X = pd.DataFrame([dict(zip(self.can_kmer, np.zeros(len(self.can_kmer))))] + [self.count_kmer(x) for x in self.dat['Seq']]).fillna(0)
        X = X.iloc[1:]
        X['Length'] = [len(x) for x in self.dat['Seq']]
        X['GC'] = [(x.count('G') + x.count('C'))/len(x) for x in self.dat['Seq']]

        y = [self.label_dict[x] for x in self.dat['Type']]

        X = X.reindex(sorted(X.columns), axis=1)

        X_train, X_test, y_train, self.y_test = train_test_split(X, y, test_size=self.test_size, random_state=self.rnd_seed)

        self.dtrain = xgb.DMatrix(X_train, label=y_train)
        self.dtest = xgb.DMatrix(X_test, label=self.y_test)

    def train(self):

        params = {
            'eta':self.eta,
            'objective':'multi:softprob',
            'eval_metric':'mlogloss',
            'nthread':self.threads,
            'num_class':len(self.incl)
        }

        print('\033[92m' + 'Using the following parameters:' + '\033[0m')
        print('eta: {}'.format(self.eta))
        print('num_rounds: {}'.format(self.num_rounds))
        print('early_stopping_rounds: {}'.format(self.early_stop))

        grid_params = [
            (max_depth, subsample, colsample_bytree)
            for max_depth in self.max_depth
            for subsample in self.subsample
            for colsample_bytree in self.colsample_bytree
        ]

        print('\033[92m' + 'Cross-validating with {} folds:'.format(self.nfold) + '\033[0m')
        
        min_mlogloss = float("Inf")
        best_params = None

        for max_depth, subsample, colsample_bytree in grid_params:
            print("CV with max_depth={}, subsample={}, colsample_bytree={}".format(
                                     max_depth,
                                     subsample,
                                     colsample_bytree))

            params['max_depth'] = max_depth
            params['subsample'] = subsample
            params['colsample_bytree'] = colsample_bytree

            cv_results = xgb.cv(
                params,
                self.dtrain,
                num_boost_round=self.num_rounds,
                seed=self.rnd_seed,
                nfold=self.nfold,
                metrics={'mlogloss'},
                early_stopping_rounds=self.early_stop
            )
            # Update best mlogloss
            mean_mlogloss = cv_results['test-mlogloss-mean'].min()
            boost_rounds = cv_results['test-mlogloss-mean'].argmin()

            print("\tmlogloss {} for {} rounds".format(mean_mlogloss, boost_rounds))

            if mean_mlogloss < min_mlogloss:
                min_mlogloss = mean_mlogloss
                best_params = (max_depth, subsample, colsample_bytree, boost_rounds)
        print("Best params: {}, {}, {}, mlogloss: {}".format(best_params[0],
                                                             best_params[1],
                                                             best_params[2],
                                                             min_mlogloss))
        

        print('\033[92m' + 'Training final model' + '\033[0m')
        
        params['max_depth'] = best_params[0]
        params['subsample'] = best_params[1]
        params['colsample_bytree'] = best_params[2]

        self.model = xgb.train(
            params,
            self.dtrain,
            num_boost_round=self.num_rounds,
            evals=[(self.dtest, "Test")],
            early_stopping_rounds=self.early_stop
        )

        self.boost_rounds = self.model.best_iteration

        self.model.save_model(self.out+'xgb_repeats.model')
        
    def test(self):

        y_pred = self.model.predict(self.dtest, iteration_range=(0, self.boost_rounds))

        conf = confusion_matrix(self.y_test, [x.argmax() for x in y_pred])
        conf_df = pd.DataFrame(conf, columns = self.incl, index = self.incl)
        conf_df.to_csv(self.out+'confusion_matrix.tab', sep='\t')

        type_acc = np.diag(conf_df) / conf_df.sum(axis=1)

        print('\033[92m' + 'Overall accuracy:' + '\033[0m')
        print(np.diag(conf).sum()/conf.sum())

        print('\033[92m' + 'Accuracy per subtype:' + '\033[0m')
        print(type_acc)

        print('\033[92m' + 'Adjusted accuracy:' + '\033[0m')
        print(type_acc.mean())

