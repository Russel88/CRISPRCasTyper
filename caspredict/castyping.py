#!/usr/bin/env python

import re
import logging
import sys

import numpy as np
import pandas as pd
import multiprocessing as mp

from itertools import chain
from scipy import ndimage

class Typer(object):

    def __init__(self, obj):
        self.master = obj
        for key, val in vars(obj).items():
            setattr(self, key, val)

    def type_operon(self, operon):
        '''
        Subtype of a single operon
        '''

        logging.debug('Typing operon '+operon)

        # Extract only the operon of interest
        tmp = self.hmm_df_all[self.hmm_df_all['operon'] == operon].sort_values('Pos')

        # If there are duplicates, keep the highest score only
        tmpX = tmp.sort_values('score', ascending=False)
        tmpX['Hmm'] = [re.sub("_.*","",x) for x in tmpX['Hmm']]
        tmpX.drop_duplicates('Hmm', inplace=True)

        start = tmp['start']
        end = tmp['end']

        start_operon = min(list(start)+list(end))
        end_operon = max(list(start)+list(end))

        # Get scores for each type
        type_scores = tmpX.iloc[:,14:].sum(axis=0)
        
        # Highest score and type with highest score
        best_score = np.amax(type_scores)
        best_type = type_scores.index.values[np.argmax(type_scores.values)]

        # At least 3 genes in operon
        if len(tmpX) >= 3:
           
            # Low score is false unless there is a signature gene
            if best_score <= 5:
                if any([x in self.signature for x in list(tmpX['Hmm'])]):
                    if sum(type_scores == best_score) > 1:
                        prediction = "Ambiguous"
                        best_type = list(type_scores.index.values[type_scores.values == np.amax(type_scores.values)])
                    else:
                        prediction = best_type
                else:
                    prediction = "False"
            
            # Solve problem with adjacent systems. Only check if at least 6 genes
            elif len(tmpX) >= 6:
                # Only types with at least one specific HMM
                zzz = tmpX.iloc[:,14:].transpose()
                zzz = zzz[zzz.apply(lambda r: any(r >= 3), axis=1)]

                # Sum of unique genes
                zzz = zzz.loc[:, zzz.apply(lambda r: sum(r > 0) == 1, axis=0)]
                zzz[zzz < 0] = 0
                zzz["sum"] = zzz.sum(axis=1)
                
                # Those with score at least 6
                zzz = zzz[zzz['sum'] >= 6]
                
                # If more than 1 type, hybrid:
                if len(zzz) > 1:
                    prediction = 'Hybrid({})'.format(','.join(zzz.index))
                # Else just choose best type
                else:
                    prediction = best_type

            else:
                # If ties, ambiguous
                if sum(type_scores == best_score) > 1:
                    prediction = "Ambiguous"
                    best_type = list(type_scores.index.values[type_scores.values == np.amax(type_scores.values)])
                # Else, choose best
                else:
                    prediction = best_type

        # 1 or 2 genes in operon
        else:

            if len(tmpX) == 2:
                first_signature = list(tmpX['Hmm'])[0] in self.signature
                second_signature = list(tmpX['Hmm'])[1] in self.signature
                accept = first_signature or second_signature
            if len(tmpX) == 1:
                accept = list(tmpX['Hmm'])[0] in self.signature
            
            # If any is a signature gene for a single gene type
            if accept:
                # If ties, ambiguous
                if sum(type_scores == best_score) > 1:
                    prediction = "Ambiguous"
                    best_type = list(type_scores.index.values[type_scores.values == np.amax(type_scores.values)])
                # If no ties
                else:
                    prediction = best_type
            # If no single gene type signature
            else:
                prediction = "False"
        

        outdict = {"Contig": list(tmp['Acc'])[0],
                   "Operon": operon,
                   "Start": start_operon,
                   "End": end_operon,
                       "Prediction": prediction,
                       "Best_type": best_type,
                       "Best_score": best_score,
                       "Genes": list(tmp['Hmm']),
                       "Positions": list(tmp['Pos']),
                       "E-values": ['{:0.2e}'.format(x) for x in list(tmp['Eval'])],
                       "CoverageSeq": [round(x,3) for x in list(tmp['Cov_seq'])],
                       "CoverageHMM": [round(x,3) for x in list(tmp['Cov_hmm'])]}
        
        return outdict


    def cluster_adj(self, data):
        '''
        Cluster adjacent genes into operons

        Params:
        data: A pandas data.frame with an Acc (accession number) column and a Pos (gene posistion) column
        dist: Int. Max allowed distance between genes in an operon

        Returns:
        List. Operon IDs with the same length and order as the input data.frame
        '''
        
        dist = self.dist
        
        positions = list(data['Pos'])
        # Create a list of zeroes to indicate positions
        pos_range = max(positions) * [0]
        # Insert ones at positions where genes are annotated
        for x in positions:
            pos_range[x-1] = 1
        # Pad
        pad = list(np.zeros(dist, dtype=int))
        pos_range_pad = pad + pos_range + pad
        # Closing to melt adjacent genes together (up to 'dist' genes between them)
        pos_range_dilated = ndimage.morphology.binary_closing(pos_range_pad, structure = list(np.ones(dist+1)))
        # Label adjacent elements
        clust_pad, nclust = ndimage.label(pos_range_dilated)
        # Remove pad
        clust = clust_pad[dist:len(clust_pad)-dist]
        # Extract cluster id for each gene
        return [list(data['Acc'])[0] + "@" + str(clust[x-1]) for x in positions]
    
    def typing(self):
        '''
        Subtyping of putative Cas operons
        '''

        if self.any_cas:
            logging.info('Subtyping putative operons')

            # Specific cut-offs
            specifics = []
            for key, value in self.cutoffs.items():
                which_sub = [i for i in list(self.hmm_df['Hmm']) if key.lower() in i.lower()]
                if len(which_sub) > 0:
                    specifics.extend(which_sub)
                    self.hmm_df = self.hmm_df[((self.hmm_df['Eval'] < float(value[0])) & 
                                                (self.hmm_df['Cov_seq'] >= float(value[1])) &
                                                (self.hmm_df['Cov_hmm'] >= float(value[2]))) |
                                                ([x not in which_sub for x in self.hmm_df['Hmm']])]

            # Apply overall thresholds for the rest
            self.hmm_df = self.hmm_df[((self.hmm_df['Cov_seq'] >= self.ocs) & 
                                        (self.hmm_df['Cov_hmm'] >= self.och) & 
                                        (self.hmm_df['Eval'] < self.oev)) |
                                        ([x in specifics for x in self.hmm_df['Hmm']])]
          
            # Define operons
            self.hmm_df = self.hmm_df.sort_values('Acc')
            operons = list(self.hmm_df.groupby('Acc').apply(self.cluster_adj))
            self.hmm_df.loc[:,'operon'] = list(chain.from_iterable(operons))

            # Load score table
            scores = pd.read_csv(self.scoring, sep=",")
            scores.fillna(0, inplace=True)

            # Signature genes for single gene types
            self.signature = [re.sub('_.*','',x) for x in list(specifics)]

            # Merge the tables
            self.hmm_df_all = pd.merge(self.hmm_df, scores, on="Hmm")

            # Assign subtype for each operon
            operons_unq = set(self.hmm_df_all['operon'])
            dictlst = [self.type_operon(operonID) for operonID in operons_unq]
            
            # Return
            self.preddf = pd.DataFrame(dictlst)

            # Check if any operons
            self.check_type()

            # Write cas operons
            self.write_type()

    def check_type(self):
        if self.any_cas:
            if len(self.preddf) == 0:
                logging.info('No operons found.')
            else:
                self.any_operon = True

    def write_type(self):
        
        if self.any_operon:
            operons_good = self.preddf[~self.preddf['Prediction'].isin(['False', 'Ambiguous'])]
            operons_put = self.preddf[self.preddf['Prediction'].isin(['False', 'Ambiguous'])]

            if len(operons_good) > 0:
                operons_good.to_csv(self.out+'cas_operons.tab', sep='\t', index=False)
            if len(operons_put) > 0:    
                operons_put.to_csv(self.out+'cas_operons_putative.tab', sep='\t', index=False)
            
