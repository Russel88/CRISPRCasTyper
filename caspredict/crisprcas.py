import logging

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
        
        # Load data
        cas = self.preddf
        cas_1 = cas[~cas['Prediction'].isin(['False', 'Ambiguous', 'Partial'])]
        crispr = pd.read_csv(self.out+'crisprs_all.tab', sep='\t')

        # Only if there is operons
        if len(cas) > 0:
            
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
                        'Subtype_Cas': list(cas_operon['Best_type'])[0],
                        'Subtype_CRISPRs': list(crispr_operon['Subtype'])
                        }

                        dicts.append(outdict)
        
        # Write
        if len(dicts) > 0:
            crispr_cas = pd.DataFrame(dicts, columns=dicts[0].keys())
            orphan_cas = cas_1[cas_1['Operon'].isin(set(cas_1['Operon']).difference(set(crispr_cas['Operon'])))]
            orphan_crispr = crispr[crispr['CRISPR'].isin(set(crispr['CRISPR']).difference(set([x for x in crispr_cas['CRISPRs'] for x in x])))]

            crispr_cas.to_csv(self.out+'CRISPR_Cas.tab', sep='\t', index=False)
            orphan_cas.to_csv(self.out+'cas_operons_orphan.tab', sep='\t', index=False)
            orphan_crispr.to_csv(self.out+'crisprs_orphan.tab', sep='\t', index=False)

