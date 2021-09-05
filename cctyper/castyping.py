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

        start = list(tmp['start'])
        end = list(tmp['end'])

        if operon in self.circ_operons:
            gene_end = np.argmax(np.diff(list(tmp['Pos'])))
            start_operon = start[gene_end+1]
            end_operon = end[gene_end]
        else:
            start_operon = min(start+end)
            end_operon = max(start+end)

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
                    best_type = list(zzz.index)
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
                # If best type is not single effector
                elif best_type not in self.single_effector:
                    prediction = "Ambiguous"
                    best_type = list(type_scores.index.values[type_scores.values > 0])
                # If no ties
                else:
                    prediction = best_type
            # If no single gene type signature
            else:
                prediction = "False"
        

        # Completion
        compl_interf = []
        compl_adapt = []
        
        # Make list of types
        if isinstance(best_type, str):
            best_type_list = [best_type]
        else:
            best_type_list = best_type
        
        for bt in best_type_list:
            # Get interference completion
            if bt in self.compl_interf.keys():
                mandatory_genes = self.compl_interf[bt]
                compl_interf.append(str(round(len([i for i in mandatory_genes if any([any([x in k for k in list(tmp['Hmm'])]) for x in i])])/len(mandatory_genes)*100))+'%')
            else:
                compl_interf.append('NA')
            # Get adaptation completion
            if bt in self.compl_adapt.keys():
                mandatory_genes = self.compl_adapt[bt]
                compl_adapt.append(str(round(len([i for i in mandatory_genes if any([any([x in k for k in list(tmp['Hmm'])]) for x in i])])/len(mandatory_genes)*100))+'%')
            else:
                compl_adapt.append('NA')

        if isinstance(best_type, str):
            compl_interf = compl_interf[0]
            compl_adapt = compl_adapt[0]

        # Strand
        interf_genes = list(chain.from_iterable(chain.from_iterable([v for k, v in self.compl_interf.items() if k in best_type_list])))
        adapt_genes = list(chain.from_iterable(chain.from_iterable([v for k, v in self.compl_adapt.items() if k in best_type_list])))

        tmp['Interf'] = tmp.apply(lambda x: any([i in x['Hmm'] for i in interf_genes]), axis=1)
        tmp['Adapt'] = tmp.apply(lambda x: any([i in x['Hmm'] for i in adapt_genes]), axis=1)

        interf_strands = set(tmp.loc[tmp['Interf'] == True, 'strand'])
        adapt_strands = set(tmp.loc[tmp['Adapt'] == True, 'strand'])

        if len(interf_strands) == 1:
            interf_strand = interf_strands.pop()
        elif len(interf_strands) == 0:
            interf_strand = "NA"
        else:
            interf_strand = 0

        if len(adapt_strands) == 1:
            adapt_strand = adapt_strands.pop()
        elif len(adapt_strands) == 0:
            adapt_strand = "NA"
        else:
            adapt_strand = 0
        
        # Collect
        outdict = {"Contig": list(tmp['Acc'])[0],
                   "Operon": operon,
                   "Start": start_operon,
                   "End": end_operon,
                       "Prediction": prediction,
                       "Complete_Interference": compl_interf,
                       "Complete_Adaptation": compl_adapt,
                       "Best_type": best_type,
                       "Best_score": best_score,
                       "Genes": list(tmp['Hmm']),
                       "Positions": list(tmp['Pos']),
                       "E-values": ['{:0.2e}'.format(x) for x in list(tmp['Eval'])],
                       "CoverageSeq": [round(x,3) for x in list(tmp['Cov_seq'])],
                       "CoverageHMM": [round(x,3) for x in list(tmp['Cov_hmm'])],
                       "Strand_Interference": interf_strand,
                       "Strand_Adaptation": adapt_strand}
        
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
        if self.circular:
            pos_range = max(self.genes[self.genes['Contig'] == list(data['Acc'])[0]]['Pos']) * [0]
        else:
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
        # Resolve circular
        is_circ = False
        if self.circular:
            if any(clust[len(clust)-dist-1:] > 0) and any(clust[:dist+1] > 0):
                last_num = clust[len(clust)-dist-1:][clust[len(clust)-dist-1:] > 0][0]
                first_num = clust[:dist+1][clust[:dist+1] > 0][0]
                if last_num != first_num:
                    clust[clust == last_num] = first_num
                    is_circ = True
        # Extract cluster id for each gene
        return [[list(data['Acc'])[0] + "@" + str(clust[x-1]) for x in positions], is_circ]
    
    def typing(self):
        '''
        Subtyping of putative Cas operons
        '''

        if self.any_cas:
            logging.info('Subtyping putative operons')

            # Save raw
            self.hmm_df_raw = self.hmm_df
            
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
            if self.circular and self.redo:
                self.genes = pd.read_csv(self.out+'genes.tab', sep='\t')
            
            self.hmm_df = self.hmm_df.sort_values('Acc')
            operons = self.hmm_df.groupby('Acc').apply(self.cluster_adj)
            self.hmm_df.loc[:,'operon'] = list(chain.from_iterable([x[0] for x in list(operons)]))

            # If any circular
            self.circ_operons = []
            any_circ = [x[1] for x in list(operons)]
            if any(any_circ):
                self.circ_operons = [sorted(i)[0] for (i, v) in zip([x[0] for x in list(operons)], any_circ) if v]

            # Prepare score table
            self.scores.fillna(0, inplace=True)
            self.cas_hmms = list(self.scores['Hmm'])
            
            # Signature genes for single gene types
            self.signature = [re.sub('_.*','',x) for x in list(specifics)]

            # Get single effector types
            single_effector_hmms = self.scores[self.scores['Hmm'].isin(list(specifics))].drop('Hmm', axis=1)
            single_effector_hmms[single_effector_hmms < 0] = 0
            self.single_effector = list(single_effector_hmms.iloc[:, single_effector_hmms.sum(axis=0).values > 0].columns)

            # Merge the tables
            self.hmm_df_all = pd.merge(self.hmm_df, self.scores, on="Hmm")

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
            
            if not self.skip_blast:
                # Get positions of Cas
                cas_positions = list(zip(self.preddf['Contig'], self.preddf['Start'], self.preddf['End']))
              
                self.flank_dict = {}
                self.flank_dict_pos = {}

                # For each contig
                for acc in set([x[0] for x in cas_positions]):
                    cas_pos = [(x[1], x[2]) for x in cas_positions if x[0] == acc]

                    # Expand by chosen distance cutoff between CRISPR and cas
                    cas_pos = [(x[0]-self.crispr_cas_dist, x[1]+self.crispr_cas_dist) for x in cas_pos] 
                    cas_pos = [x if x[0]>0 else (1, x[1]) for x in cas_pos]

                    # Check if any overlap, and merge if true
                    if len(cas_pos) > 1:
                        def recursive_merge(ll, start_index = 0):
                            for i in range(start_index, len(ll) - 1):
                                if ll[i][1] > ll[i+1][0]:
                                    new_start = ll[i][0]
                                    new_end = ll[i+1][1]
                                    ll[i] = (new_start, new_end)
                                    del ll[i+1]
                                    return recursive_merge(ll.copy(), start_index=i)
                            return ll

                        cas_pos_sort = sorted(cas_pos)
                        cas_pos = recursive_merge(cas_pos_sort.copy())

                    n = 0
                    for i in cas_pos:
                        n += 1
                        with open(self.out+'Flank.fna', 'a') as handle:
                            handle.write('>'+acc+'-'+str(n)+'\n')
                            handle.write(str(self.seq_dict[acc][(i[0]-1):i[1]])+'\n')
                            self.flank_dict[acc+'-'+str(n)] = str(self.seq_dict[acc][(i[0]-1):i[1]])
                            self.flank_dict_pos[acc+'-'+str(n)] = (i[0], i[1])
                            
