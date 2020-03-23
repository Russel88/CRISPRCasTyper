import logging
import re

import pandas as pd

class CRISPRCas(object):
    
    def __init__(self, obj):
        self.master = obj
        for key, val in vars(obj).items():
            setattr(self, key, val)

    def crisprcas(self):
        
        # Define distance functions
        def dist(x,y):
            return y[0]-x[1] if y[0]>x[1] else x[0]-y[1]

        def dist_ll(x,ll):
            return [dist(x,y) for y in ll]
        
        if not self.any_crispr:
            self.crisprsall = []
        
        # Only if there is operons and crisprs
        if self.any_operon and self.any_crispr:

            logging.info('Connecting Cas operons and CRISPR arrays')

            # Load data
            cas = self.preddf
            cas = cas[~cas['Prediction'].isin(['False'])]
            cas_1 = cas[~cas['Prediction'].isin(['False', 'Ambiguous', 'Partial'])]
            crispr = pd.read_csv(self.out+'crisprs_all.tab', sep='\t')
            self.crisprsall = crispr

            dicts = []
            # Loop over contig
            for contig in set(cas['Contig']):
                cas_sub = cas[cas['Contig'] == contig]
                crispr_sub = crispr[crispr['Contig'] == contig]
                
                # Loop over operons
                for operon in set(cas_sub['Operon']):
                    cas_operon = cas_sub[cas_sub['Operon'] == operon]
                    
                    # Find distances between operon and crisprs
                    dists = dist_ll((int(cas_operon['Start']), int(cas_operon['End'])), zip(crispr_sub['Start'], crispr_sub['End']))
                    
                    # Only crisprs closer than dist threshold
                    crispr_operon = crispr_sub[[x <= self.crispr_cas_dist for x in dists]]
                    
                    if len(crispr_operon) > 0:
                        outdict = {
                        'Contig': list(cas_operon['Contig'])[0],
                        'Operon': list(cas_operon['Operon'])[0],
                        'Operon_Pos': [list(cas_operon['Start'])[0], list(cas_operon['End'])[0]],     
                        'CRISPRs': list(crispr_operon['CRISPR']),
                        'Distances': [0 if x<0 else x for x in dists if x <= self.crispr_cas_dist],
                        'Prediction_Cas': list(cas_operon['Prediction'])[0],
                        'Prediction_CRISPRs': list(crispr_operon['Prediction']),
                        'Subtype_Cas': list(cas_operon['Best_type'])[0],
                        'Subtype_CRISPRs': list(crispr_operon['Subtype'])
                        }

                        dicts.append(outdict)

            # Consensus 
            if len(dicts) > 0:
                crispr_cas = pd.DataFrame(dicts, columns=dicts[0].keys())
                self.orphan_cas = cas_1[cas_1['Operon'].isin(set(cas_1['Operon']).difference(set(crispr_cas['Operon'])))]
                self.orphan_crispr = crispr[crispr['CRISPR'].isin(set(crispr['CRISPR']).difference(set([x for x in crispr_cas['CRISPRs'] for x in x])))]
            
                pred_lst = []
                for index, row in crispr_cas.iterrows():
                    # Choose nearest CRISPR prediction
                    Prediction_CRISPR = row['Prediction_CRISPRs'][row['Distances'].index(min(row['Distances']))]
                    
                    # Get Cas predictions
                    Prediction_Cas = row['Prediction_Cas']
                    Best_Cas = row['Subtype_Cas']

                    # If agree
                    if Prediction_Cas == Prediction_CRISPR:
                        Prediction = Prediction_Cas
                    # If not agree
                    else:
                        # If Cas ambiguous
                        if Prediction_Cas == "Ambiguous":
                            if Prediction_CRISPR in Best_Cas:
                                Prediction = Prediction_CRISPR
                            # If overall type match
                            elif re.sub('-.*$', '', Prediction_CRISPR) in [re.sub('-.*$', '', x) for x in Best_Cas]:
                                Prediction = re.sub('-.*$', '', Prediction_CRISPR)
                            else:
                                Prediction = "Unknown"
                        # If Cas not False or Partial
                        elif Prediction_Cas not in ('False', 'Partial'):
                            Prediction = Prediction_Cas
                        # If Cas False or Partial    
                        else:
                            if Best_Cas == Prediction_CRISPR:
                                Prediction = Best_Cas+'(Partial)'
                            # If overall type match
                            elif re.sub('-.*$', '', Best_Cas) == re.sub('-.*$', '', Prediction_CRISPR):
                                Prediction = re.sub('-.*$', '', Prediction_CRISPR)+'(Partial)'
                            else:
                                Prediction = "Unknown"
                                
                    pred_lst.append(Prediction)

                # Insert predictions
                crispr_cas['Prediction'] = pred_lst

                # Subset columns
                self.crispr_cas = crispr_cas[['Contig', 'Operon', 'Operon_Pos', 'Prediction', 'CRISPRs', 'Distances', 'Prediction_Cas', 'Prediction_CRISPRs']]
            
            # Write
            if len(dicts) > 0:

                # Split Crispr cas in putative and good
                crispr_cas_good = self.crispr_cas[~self.crispr_cas['Prediction'].str.contains('Unknown|Partial')]
                crispr_cas_put = self.crispr_cas[self.crispr_cas['Prediction'].str.contains('Unknown|Partial')]

                if len(crispr_cas_good) > 0:
                    crispr_cas_good.to_csv(self.out+'CRISPR_Cas.tab', sep='\t', index=False)
                if len(crispr_cas_put) > 0:
                    crispr_cas_put.to_csv(self.out+'CRISPR_Cas_putative.tab', sep='\t', index=False)
                if len(self.orphan_cas) > 0:
                    self.orphan_cas.to_csv(self.out+'cas_operons_orphan.tab', sep='\t', index=False)
                if len(self.orphan_crispr) > 0:
                    self.orphan_crispr.to_csv(self.out+'crisprs_orphan.tab', sep='\t', index=False)
                

