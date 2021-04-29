import os
import subprocess
import logging
import sys
import re
import math
import random

import pandas as pd
import statistics as st

from Bio import pairwise2
from Bio import SeqIO
from Bio.Seq import Seq
from joblib import Parallel, delayed

from cctyper.minced import CRISPR

class RepeatMatch(object):
    
    def __init__(self, obj):
        self.master = obj
        for key, val in vars(obj).items():
            setattr(self, key, val)
    
    def run(self):
        
        if self.any_operon and not self.skip_blast:
            
            logging.info("BLASTing for CRISPR near cas operons")
            
            self.make_db()
            self.align()
            self.clust()

    def make_db(self):
        '''
        Make a BLAST database
        '''

        logging.debug('Making BLAST database')

        subprocess.run(['makeblastdb', 
                        '-dbtype', 'nucl', 
                        '-in', self.out+'Flank.fna',
                        '-out', self.out+'Flank'], 
                        stdout=subprocess.DEVNULL)
    
    def align(self):
        '''
        BLASTing repeat database
        '''

        logging.debug('BLASTing repeats')

        # BLASTn
        subprocess.run(['blastn', 
                        '-task', 'blastn-short', 
                        '-query', self.repeatdb,
                        '-db', self.out+'Flank',
                        '-outfmt', '6',
                        '-out', self.out+'blast.tab',
                        '-num_threads', str(self.threads),
                        '-perc_identity', str(90),
                        '-qcov_hsp_perc', str(90)])


    def clust(self):
        '''
        Load the BLAST table and run the different clustering steps
        '''

        logging.debug('Clustering matches into arrays')
        
        # Load blast table and add lengths
        self.df = pd.read_csv(self.out+'blast.tab', sep='\t', header=None,
            names=('Repeat', 'Acc', 'Identity', 'Alignment', 'Mismatches', 'Gaps',
                   'Repeat_start', 'Repeat_end', 'Acc_start', 'Acc_end', 'Evalue', 'Score'))
       
        # If any matches
        if len(self.df) > 0:
            # Create new columns
            self.df['Min'] = [min(x,y) for x,y in zip(self.df['Acc_start'],self.df['Acc_end'])]
            self.df['Max'] = [max(x,y) for x,y in zip(self.df['Acc_start'],self.df['Acc_end'])]

            # Keep only best matches if overlapping
            self.remove_overlap()

            # Add repeats
            self.add_repeats()

            # Cluster matches in arrays
            self.cluster_adj()

            count_dict = self.df_cluster.groupby('Cluster')['Cluster'].count().to_dict()
            
            # At least 3 repeats
            cluster_array = [x for x in count_dict if count_dict[x] > 2]
         
            # If any arrays
            if len(cluster_array) > 0:

                self.df_array = self.df_cluster[self.df_cluster['Cluster'].isin(cluster_array)]
                self.convert_array()

    def overlap(self,x,y):
        '''
        Evaluate whether two (start,end) tuples overlap
        '''
        return x[0] <= y[1] and y[0] <= x[1]

    def overlap_any(self,x,ll):
        '''
        Evaluate whether a (start,end) tuple has overlap with any in a list of tuples
        '''
        return any([self.overlap(x,y) for y in ll])

    def distance(self,x,y):
        '''
        Calculate distance between two (start,end) tuples
        '''
        return y[0]-x[1] if y[0]>x[1] else x[0]-y[1]

    def distance_all(self,x,ll):
        '''
        Get all distances between a (start,end) tuple and a list of tuples
        '''
        return [self.distance(x,y) for y in ll]
    
    def remove_overlap(self):
        '''
        If matches overlap keep only the best
        '''

        logging.debug('Removing overlapping matches')

        # Sort by alignment quality
        self.df = self.df.sort_values(['Acc', 'Score'], ascending=False) 

        overlap_lst = []
        for i in set(self.df['Acc']):
            tmp = self.df[self.df['Acc'] == i]
            
            # First remove those with similar start or end
            tmp = tmp.drop_duplicates('Min')
            tmp = tmp.drop_duplicates('Max')
            tmp = tmp.drop_duplicates('Acc_start')
            tmp = tmp.drop_duplicates('Acc_end')

            # Then traverse through matches comparing only with previous
            pos = tmp[['Min','Max']].values
            keep = []
            matches_all = []
            # For each match
            for ind, k in enumerate(pos):
                # If no overlaps with any previous, keep
                if not self.overlap_any(k, matches_all):
                    keep.append(ind)
                    matches_all.append(k)

            overlap_lst.append(tmp.iloc[keep,:])
        
        # If several contigs, concatenate
        self.df_overlap = pd.concat(overlap_lst)

    def cluster_adj(self):
        '''
        Cluster adjacent matches into arrays
        '''

        logging.debug('Clustering matches')

        # Sort by position
        self.df_overlap = self.df_overlap.sort_values('Min')

        cluster_df_lst = []
        cluster = 0
        # For each contig
        for i in set(self.df_overlap['Acc']):
            tmp = self.df_overlap[self.df_overlap['Acc'] == i]

            pos = tmp[['Min','Max']].values
            cluster_list = []
            # Loop over complete matches
            for ind, k in enumerate(pos):
                # Keep first match
                if ind == 0:
                    cluster_list.append(cluster)
                    arrays_cluster = [k]
                else:
                    # If match within Xbp of any previous, add match to current cluster
                    if min(self.distance_all(k, arrays_cluster)) <= 100:
                        cluster_list.append(cluster)
                        arrays_cluster.append(k)
                    # If match > Xbp from previous, initiate new cluster
                    else:
                        cluster += 1
                        arrays_cluster = [k]
                        cluster_list.append(cluster)
            
            tmp.insert(len(tmp.columns), 'Cluster', cluster_list)
            cluster_df_lst.append(tmp)
            
            # Increment cluster ID for next acc
            cluster += 1

        # If several contigs, concatenate
        self.df_cluster = pd.concat(cluster_df_lst)

    def get_sequence(self, acc, start, end):
        '''
        Return sequence from position information
        '''

        if start < 1:
            start = 1
        
        return(''.join(self.flank_dict[acc][(start-1):end]))

    def add_repeats(self):
        '''
        Add repeats to the no-overlap dataframe
        '''
        
        logging.debug('Adding repeat sequences')

        self.df_overlap.insert(len(self.df_overlap.columns), 'Sequence', self.df_overlap.apply(lambda row: self.get_sequence(row['Acc'], row['Min'], row['Max']), axis=1))

    def convert_array(self):
        '''
        Convert array to CRISPR object
        '''

        logging.debug('Converting array dataframe')

        crisprs = []
        # For each array
        cls = set(self.df_array['Cluster'])
        for cl in cls:
            tmp = self.df_array[self.df_array['Cluster'] == cl]
            acc = list(tmp['Acc'])[0]

            # Initiate CRISPR object
            crisp_tmp = CRISPR(acc, self.exact_stats)

            # Get positions
            flank_start = self.flank_dict_pos[acc][0]-1
            crisp_tmp.setPos(str(min(tmp['Min'])+flank_start),  str(max(tmp['Max'])+flank_start))
           
            # Add repeats
            for repseq in list(tmp['Sequence']):
                crisp_tmp.addRepeat(repseq)

            # Get spacers
            spacers = [str(self.get_sequence(acc, x[0]+1, x[1]-1)) for x in zip(tmp['Max'][:(len(tmp)-1)], tmp['Min'][1:])]
            for spaseq in spacers:
                crisp_tmp.addSpacer(spaseq)

            # Stats
            crisp_tmp.getConsensus()
            crisp_tmp.stats(self.threads, self.repeat_id, self.spacer_id, self.spacer_sem)
            
            # Edit naming
            crisp_tmp.sequence = re.sub('-[0-9]*$', '', crisp_tmp.sequence)
            
            # Save the instance
            crisprs.append(crisp_tmp)

        # Remove very irregular arrays
        crisprs_good = []
        for crisp_tmp in crisprs:
            ok = 0
            if crisp_tmp.spacer_identity < 60:
                ok += 1
            if crisp_tmp.repeat_identity > 65:
                ok += 1
            if crisp_tmp.spacer_sem < 6:
                ok += 1
            if ok == 3:
                crisprs_good.append(crisp_tmp)
        crisprs = crisprs_good
        
        # Continue if any left
        if len(crisprs) > 0:

            # Update CRISPR object
            crisprs_new = []
            if hasattr(self, 'crisprs'):
                for crisp in crisprs:
                    # Check if overlapping
                    crisp_in_seq = [x for x in self.crisprs if x.sequence == crisp.sequence]
                    if not any([x.start <= crisp.end and x.end >= crisp.start for x in crisp_in_seq]):
                        # Append if not overlapping
                        self.crisprs.append(crisp)
                        crisprs_new.append(crisp)
            else:
                self.crisprs = crisprs
            crisprs = crisprs_new
            
            header = '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format('Contig',
                                                                        'CRISPR',
                                                                        'Start',
                                                                        'End',
                                                                        'Consensus_repeat',
                                                                        'N_repeats',
                                                                        'Repeat_len',
                                                                        'Spacer_len_avg',
                                                                        'Repeat_identity',
                                                                        'Spacer_identity',
                                                                        'Spacer_len_sem',
                                                                        'Trusted')
           
            def write_crisp(handle, cris):
                handle.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(cris.sequence,
                                                       cris.crispr,
                                                       cris.start,
                                                       cris.end,
                                                       cris.cons,
                                                       len(cris.repeats),
                                                       cris.repeat_len,
                                                       cris.spacer_len,
                                                       cris.repeat_identity,
                                                       cris.spacer_identity,
                                                       cris.spacer_sem,
                                                       cris.trusted))
        
            # Update CRISPR file
            if os.path.exists(self.out+'crisprs_all.tab'):
                f = open(self.out+'crisprs_all.tab', 'a')
            else:
                f = open(self.out+'crisprs_all.tab', 'w')
                f.write(header)
            for crisp in crisprs:
                write_crisp(f, crisp)
            f.close()

            # Write spacers
            if not os.path.exists(self.out+'spacers'):
                os.mkdir(self.out+'spacers')

            for crisp in crisprs:
                f = open(self.out+'spacers/{}.fa'.format(crisp.crispr), 'w')
                n = 0
                for sq in crisp.spacers:
                    n += 1
                    f.write('>{}:{}\n'.format(crisp.crispr, n))
                    f.write('{}\n'.format(sq))
                f.close()

