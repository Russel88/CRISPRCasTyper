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
           
            # Solve problem with adjacent systems. Only check if at least 6 genes
            if len(tmpX) >= 6:
                # Only types with at least one specific HMM
                zzz = tmpX.iloc[:,14:].transpose()
                zzz= zzz[zzz.apply(lambda r: any(r >= 3), axis=1)]

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


        # 2 genes in operon
        elif len(tmpX) == 2:

            # Both somewhat good matches and score higher than 4
            first_gene_good = float(list(tmpX['Cov_seq'])[0]) >= self.tcs and float(list(tmpX['Cov_hmm'])[0]) >= self.tch and float(list(tmpX['Eval'])[0]) < self.tev
            second_gene_good = float(list(tmpX['Cov_seq'])[1]) >= self.tcs and float(list(tmpX['Cov_hmm'])[0]) >= self.tch and float(list(tmpX['Eval'])[1]) < self.tev

            if (first_gene_good and second_gene_good) and best_score >= 4:
                # If ties, ambiguous
                if sum(type_scores == best_score) > 1:
                    prediction = "Ambiguous"
                    best_type = "Ambiguous"
                # If no ties
                else:
                    # If one gene subtype, then assign
                    if best_type in self.single_gene_types:
                        prediction = best_type
                    # Not one-gene subtype system hit = Partial
                    else:
                        prediction = "Partial"
            # Low score or low quality double gene operon = Trash
            else:
                prediction = "False"


        # Only 1 gene
        else:

            # Only high quality
            lonely_gene_good = float(tmpX['Cov_seq']) >= self.scs and float(tmpX['Cov_hmm']) >= self.sch and float(tmpX['Eval']) < self.sev

            if lonely_gene_good and best_score >= 4:
                # If ties, ambiguous
                if sum(type_scores == best_score) > 1:
                    prediction = "Ambiguous"
                    best_type = "Ambiguous"
                # If no ties
                else:
                    # One gene subtype
                    if best_type in self.single_gene_types:
                        prediction = best_type
                    # High quality, lonely, not one-gene subtype system hit = Partial
                    else:
                        prediction = "Partial"
            # Low quality, single gene operon = Trash
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


    def typing(self):
        '''
        Subtyping of putative Cas operons
        '''

        if self.any_cas:
            logging.info('Subtyping putative operons')
            
            # Apply overall thresholds
            self.hmm_df = self.hmm_df[(self.hmm_df['Cov_seq'] >= self.ocs) & 
                                        (self.hmm_df['Cov_hmm'] >= self.och) & 
                                        (self.hmm_df['Eval'] < self.oev)]

            # V-F specific thresholds
            VF = [i for i in list(self.hmm_df['Hmm']) if 'cas12f' in i or 'cas12_x' in i]
            
            self.hmm_df = self.hmm_df[((self.hmm_df['Cov_hmm'] >= self.vfc) & 
                                        (self.hmm_df['Eval'] < self.vfe)) |
                                        ([x not in VF for x in self.hmm_df['Hmm']])]


            # Define operons
            # First define function for finding them
            def cluster_adj(data, dist=self.dist):
                '''
                Cluster adjacent genes into operons

                Params:
                data: A pandas data.frame with an Acc (accession number) column and a Pos (gene posistion) column
                dist: Int. Max allowed distance between genes in an operon

                Returns:
                List. Operon IDs with the same length and order as the input data.frame
                '''
                
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
            
            self.hmm_df.sort_values('Acc', inplace=True)
            operons = list(self.hmm_df.groupby('Acc').apply(cluster_adj))
            self.hmm_df['operon'] = list(chain.from_iterable(operons))

            # Load score table
            scores = pd.read_csv(self.scoring, sep=",")
            scores.fillna(0, inplace=True)

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
            operons_good = self.preddf[~self.preddf['Prediction'].isin(['False', 'Ambiguous', 'Partial'])]
            operons_put = self.preddf[self.preddf['Prediction'].isin(['False', 'Ambiguous', 'Partial'])]

            if len(operons_good) > 0:
                operons_good.to_csv(self.out+'cas_operons.tab', sep='\t', index=False)
            if len(operons_put) > 0:    
                operons_put.to_csv(self.out+'cas_operons_putative.tab', sep='\t', index=False)
            
