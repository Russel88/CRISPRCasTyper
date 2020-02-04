#!/usr/bin/env python

import multiprocessing as mp
from itertools import chain
import numpy as np
import pandas as pd
from scipy import ndimage
import re

def cluster_adj(data, dist):
    '''
    Cluster adjacent genes into operons
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


def type_operon(dat_all, operon, one_gene_eval, one_gene_cov_seq, one_gene_cov_hmm, two_gene_eval, two_gene_cov_seq, two_gene_cov_hmm, single_gene_types):
    '''
    Subtype of a single operon
    '''

    tmp = dat_all[dat_all['operon'] == operon].sort_values('Pos')

    # If there are duplicates, keep the highest score only
    tmpX = tmp.sort_values('score', ascending=False)
    tmpX['Hmm'] = [re.sub("_.*","",x) for x in tmpX['Hmm']]
    tmpX.drop_duplicates('Hmm', inplace=True)

    start = tmp['start']
    end = tmp['end']

    start_operon = min(list(start)+list(end))
    end_operon = max(list(start)+list(end))

    type_scores = tmpX.iloc[:,13:].sum(axis=0)

    best_score = np.amax(type_scores)
    best_type = type_scores.index.values[np.argmax(type_scores.values)]


    # At least 3 genes in operon
    if len(tmpX) >= 3:
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
        first_gene_good = float(list(tmpX['Cov_seq'])[0]) >= two_gene_cov_seq and float(list(tmpX['Cov_hmm'])[0]) >= two_gene_cov_hmm and float(list(tmpX['Eval'])[0]) < two_gene_eval
        second_gene_good = float(list(tmpX['Cov_seq'])[1]) >= two_gene_cov_seq and float(list(tmpX['Cov_hmm'])[0]) >= two_gene_cov_hmm and float(list(tmpX['Eval'])[1]) < two_gene_eval

        if (first_gene_good and second_gene_good) and best_score >= 4:
            # If ties, ambiguous
            if sum(type_scores == best_score) > 1:
                prediction = "Ambiguous"
                best_type = "Ambiguous"
            # If no ties
            else:
                # If one gene subtype, then assign
                if best_type in single_gene_types:
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
        lonely_gene_good = float(tmpX['Cov_seq']) >= one_gene_cov_seq and float(tmpX['Cov_hmm']) >= one_gene_cov_hmm and float(tmpX['Eval']) < one_gene_eval

        if lonely_gene_good and best_score >= 4:
            # If ties, ambiguous
            if sum(type_scores == best_score) > 1:
                prediction = "Ambiguous"
                best_type = "Ambiguous"
            # If no ties
            else:
                # One gene subtype
                if best_type in single_gene_types:
                    prediction = best_type
                # High quality, lonely, not one-gene subtype system hit = Partial
                else:
                    prediction = "Partial"
        # Low quality, single gene operon = Trash
        else:
            prediction = "False"


    outdict = {"Operon": operon,
               'Start': start_operon,
               'End': end_operon,
                   "Prediction": prediction,
                   "Best_type": best_type,
                   "Best_score": best_score,
                   "Genes": list(tmp['Hmm']),
                   "Positions": list(tmp['Pos']),
                   "E-values": ['{:0.2e}'.format(x) for x in list(tmp['Eval'])],
                   "CoverageSeq": [round(x,3) for x in list(tmp['Cov_seq'])],
                   "CoverageHMM": [round(x,3) for x in list(tmp['Cov_hmm'])]}
    
    return outdict


def typing(dat, out, threads, scoring, dist, overall_eval, overall_cov_seq, overall_cov_hmm, one_gene_eval, one_gene_cov_seq, one_gene_cov_hmm, two_gene_eval, two_gene_cov_seq, two_gene_cov_hmm, single_gene_types, VF_eval, VF_cov_hmm):
    '''
    Subtyping of putative Cas operons
    '''

    # Apply overall thresholds
    dat = dat[(dat['Cov_seq'] >= overall_cov_seq) & (dat['Cov_hmm'] >= overall_cov_hmm) & (dat['Eval'] < overall_eval)]

    # V-F specific thresholds
    VF = [i for i in list(dat['Hmm']) if 'cas12f' in i or 'cas12_x' in i]
    
    dat = dat[((dat['Cov_hmm'] >= VF_cov_hmm) & (dat['Eval'] < VF_eval)) |
              ([x not in VF for x in dat['Hmm']])]


    # Define operons
    dat.sort_values('Acc', inplace=True)
    dat['operon'] = list(chain.from_iterable(list(dat.groupby('Acc').apply(cluster_adj, dist))))

    # Load score table
    scores = pd.read_csv(scoring, sep=",")
    scores.fillna(0, inplace=True)

    # Merge the tables
    dat_all = pd.merge(dat, scores, on="Hmm")

    # Start multiprocess
    pool = mp.Pool(threads)

    # Assign subtype for each operon
    dictlst = [pool.apply(type_operon, args=(dat_all, operonID, one_gene_eval, one_gene_cov_seq, one_gene_cov_hmm, two_gene_eval, two_gene_cov_seq, two_gene_cov_hmm, single_gene_types)) for operonID in set(dat_all['operon'])]

    # Close multiprocess
    pool.close()

    # Return
    preddf = pd.DataFrame(dictlst)
    return(preddf)


