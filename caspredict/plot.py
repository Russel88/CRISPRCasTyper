import os
import logging
import re

import pandas as pd

from PIL import Image
from PIL import ImageDraw
from PIL import ImageFont

class Map(object):
    
    def __init__(self, obj):
        self.master = obj
        for key, val in vars(obj).items():
            setattr(self, key, val)
        
        self.font = ImageFont.truetype(os.path.join(self.db, 'arial.ttf'), 30)
        self.fontB = ImageFont.truetype(os.path.join(self.db, 'arial.ttf'), 50)
        self.fontS = ImageFont.truetype(os.path.join(self.db, 'arial.ttf'), 20)

    def draw_gene(self, start, end, strand, name, n, z, col):
        if strand > 0:
            self.draw.polygon(((self.scale/50*start, n*20*self.scale),
                          (self.scale/50*end-5*self.scale, n*20*self.scale),
                          (self.scale/50*end, 2.5*self.scale+n*20*self.scale),
                          (self.scale/50*end-5*self.scale, 5*self.scale+n*20*self.scale),
                          (self.scale/50*start, 5*self.scale+n*20*self.scale)),
                     fill=col, outline=(255, 255, 255))
        else:
            self.draw.polygon(((self.scale/50*start+5*self.scale, n*20*self.scale),
                          (self.scale/50*end, n*20*self.scale),
                          (self.scale/50*end, 5*self.scale+n*20*self.scale),
                          (self.scale/50*start+5*self.scale, 5*self.scale+n*20*self.scale),
                          (self.scale/50*start, 2.5*self.scale+n*20*self.scale)),
                     fill=col, outline=(255, 255, 255))
        
        if isinstance(name, str):
            name = re.sub('_[0-9]*_.*', '', name)
        else:
            name = str(name)

        if z % 2 == 1:
            self.draw.text((self.scale/50*start+5, n*20*self.scale-3*self.scale), name.title(), (0,0,0), font=self.font)
        else:
            self.draw.text((self.scale/50*start+5, n*20*self.scale+5*self.scale), name.title(), (0,0,0), font=self.font)

    def draw_array(self, start, end, subtype, n, z):
        self.draw.polygon(((self.scale/50*start, n*20*self.scale), 
                      (self.scale/50*end, n*20*self.scale), 
                      (self.scale/50*end, 5*self.scale+n*20*self.scale), 
                      (self.scale/50*start, 5*self.scale+n*20*self.scale)), 
                     fill=(0, 0, 255), outline=(255, 255, 255))
        self.draw.text((self.scale/50*start+self.scale/10, n*20*self.scale-3*self.scale), subtype, (0,0,0), font=self.font)

    def draw_name(self, n, pred, contig, start, end):
        self.draw.text((self.scale/10, n*20*self.scale-10*self.scale), '{}: {}({}-{})'.format(pred, contig, start, end),
                  (0,0,0), font=self.fontB)

    def draw_system(self, cas, crispr, n):
        z = 0
        if len(cas) > 0:
            for i in cas:
                z += 1
                self.draw_gene(i[0], i[1], i[2], i[3], n, z, i[4])
        if len(crispr) > 0:
            for i in crispr:
                z += 1
                self.draw_array(i[0], i[1], i[2], n, z)

    def criscas_len(self, cc, cca):
        ccs = [cca[cca['CRISPR'].isin(x)] for x in cc['CRISPRs']]
        startsO = [x[0] for x in cc['Operon_Pos']]
        endsO = [x[1] for x in cc['Operon_Pos']]
        startsC = [min(x['Start']) for x in ccs]
        endsC = [max(x['End']) for x in ccs]
        return [x[0]-x[1] for x in zip([max(x) for x in zip(endsO, endsC)], [min(x) for x in zip(startsO, startsC)])]

    def get_longest(self, crisO, casO, criscas, crisA):
        crisO_M, casO_M, cc1_M = 0, 0, 0
        if len(crisO) > 0: crisO_M = max(crisO['End']-crisO['Start'])
        if len(casO) > 0: casO_M = max(casO['End']-casO['Start'])
        if len(criscas) > 0: cc1_M = max(self.criscas_len(criscas, crisA))
        return(max(crisO_M, casO_M, cc1_M))

    def expandCas(self, contig, pos, startPos, array=False):
        first_cas = min(pos)
        last_cas = max(pos)
        if array:
            missing_cas = list(range(first_cas-self.expand+1, last_cas+self.expand))
        else:
            incl_cas = list(range(first_cas-self.expand, last_cas+self.expand+1))
            missing_cas = [x for x in incl_cas if x not in pos]
        
        add_these = self.genes[(self.genes['Contig'] == contig) & (self.genes['Pos'].isin(missing_cas))]
        
        add_starts = [self.expand*self.plotexpand + 1 + x - startPos for x in list(add_these['Start'])]
        add_ends = [self.expand*self.plotexpand + 1 + x - startPos for x in list(add_these['End'])]

        names = list(add_these['Pos'])
        cols = list(((150,150,150),)*len(add_starts))

        # Add putative
        add_putative = [x in list(self.hmm_df['Pos']) for x in list(add_these['Pos'])]
        casNames = [list(self.hmm_df[self.hmm_df['Pos'] == x]['Hmm']) for x in list(add_these['Pos'])]
        casNames = [x[0] if len(x)>0 else x for x in casNames]
        cols = [x[0] if not x[1] else (0,150,0) for x in zip(cols, add_putative)]
        names = [x[0] if not x[1] else x[2] for x in zip(names, add_putative, casNames)]

        expand_list = list(zip(add_starts,
                          add_ends,
                          list(add_these['Strand']),
                          names,
                          cols))
        
        return expand_list    

    def plot(self):

        total = 0

        # Combine orphan and ambiguous cas operons
        if self.any_operon:
            cas_ambi = self.preddf[self.preddf['Prediction'] == 'Ambiguous']
            try:
                casAmbiOrph = pd.concat([self.orphan_cas, cas_ambi])
            except:
                cas_good = self.preddf[~self.preddf['Prediction'].isin(['False','Ambiguous','Partial'])]
                casAmbiOrph = pd.concat([cas_good, cas_ambi])
            total += len(casAmbiOrph)

        try:
            total += len(self.orphan_crispr)
        except:
            self.orphan_crispr = self.crisprsall
            total += len(self.orphan_crispr)

        try:
            total += len(self.crispr_cas)
        except:
            self.crispr_cas = []

        # Plot
        if (not self.noplot) and total > 0:
            
            logging.info('Plotting map of CRISPR-Cas loci')
            
            width = self.get_longest(self.orphan_crispr, casAmbiOrph, self.crispr_cas, self.crisprsall) 

            self.genes = pd.read_csv(self.out+'genes.tab', sep='\t') 
            width = width + (self.plotexpand * self.expand * 2)

            self.im = Image.new('RGB', (int(round(self.scale/50*width+self.scale*10)), int(round((total+1)*20*self.scale))), (255, 255, 255))
            self.draw = ImageDraw.Draw(self.im)

            if not self.nogrid:
                # Draw grid
                y_start = 8*self.scale
                y_end = self.im.height
                step_size = int(round(1000*self.scale/50))

                for x in range(step_size, self.im.width, step_size):
                    line = ((x, y_start), (x, y_end))
                    self.draw.line(line, fill=(150,150,150), width=int(self.scale/20))
                    self.draw.text((x-self.scale*4, self.scale*5), str(int(x/(self.scale/50))), (100,100,100), font=self.fontS)

            k = 0
            # Draw CRISPR-Cas
            if len(self.crispr_cas) > 0:
                for i in list(self.crispr_cas['Operon']):
                    k += 1
                    logging.debug('Plotting '+i)

                    # Get data
                    contig = list(self.crispr_cas[self.crispr_cas['Operon'] == i]['Contig'])[0]
                    prediction = list(self.crispr_cas[self.crispr_cas['Operon'] == i]['Prediction'])[0]

                    # Cas
                    posCas = list(self.preddf[self.preddf['Operon'] == i]['Positions'])[0]
                    nameCas = list(self.preddf[self.preddf['Operon'] == i]['Genes'])[0]
                    hmmSub = self.hmm_df[self.hmm_df['Acc'] == contig]
                    startsCas = [list(hmmSub[hmmSub['Pos'] == x]['start'])[0] for x in posCas]
                    endsCas = [list(hmmSub[hmmSub['Pos'] == x]['end'])[0] for x in posCas]
                    strands = [list(hmmSub[hmmSub['Pos'] == x]['strand'])[0] for x in posCas]

                    # Crisprs
                    crisprs = list(self.crispr_cas[self.crispr_cas['Operon'] == i]['CRISPRs'])[0]
                    startsCris = [list(self.crisprsall[self.crisprsall['CRISPR'] == x]['Start'])[0] for x in crisprs]
                    endsCris = [list(self.crisprsall[self.crisprsall['CRISPR'] == x]['End'])[0] for x in crisprs]
                    nameCris = [list(self.crisprsall[self.crisprsall['CRISPR'] == x]['Prediction'])[0] for x in crisprs]

                    # Find start of loci
                    startPos = min(startsCas + startsCris)

                    # Draw name
                    self.draw_name(k, prediction, i, startPos, max(endsCas + endsCris))
                    
                    # Adjust positions
                    startsCas = [self.expand*self.plotexpand + 1 + x - startPos for x in startsCas]
                    endsCas = [self.expand*self.plotexpand + 1 + x - startPos for x in endsCas]
                    startsCris = [self.expand*self.plotexpand + 1 + x - startPos for x in startsCris]
                    endsCris = [self.expand*self.plotexpand + 1 + x - startPos for x in endsCris]
                    
                    # Draw
                    cas_list = list(zip(startsCas, endsCas, strands, nameCas, list(((255,0,0),)*len(nameCas))))

                    # Expand
                    expand_list = self.expandCas(contig, posCas, startPos)
                    cas_list = cas_list + expand_list

                    cas_list = sorted(cas_list, key=lambda x: x[0])
                    self.draw_system(cas_list, list(zip(startsCris, endsCris, nameCris)), k)
            
            # Draw Orphan and Ambibguous Cas
            if len(casAmbiOrph) > 0:
                for i in list(casAmbiOrph['Operon']):
                    k += 1
                    logging.debug('Plotting '+i)
                    
                    # Get data
                    contig = list(casAmbiOrph[casAmbiOrph['Operon'] == i]['Contig'])[0]
                    pos = list(casAmbiOrph[casAmbiOrph['Operon'] == i]['Positions'])[0]
                    casName = list(casAmbiOrph[casAmbiOrph['Operon'] == i]['Genes'])[0]
                    hmmSub = self.hmm_df[self.hmm_df['Acc'] == contig]
                    starts = [list(hmmSub[hmmSub['Pos'] == x]['start'])[0] for x in pos]
                    ends = [list(hmmSub[hmmSub['Pos'] == x]['end'])[0] for x in pos]
                    strands = [list(hmmSub[hmmSub['Pos'] == x]['strand'])[0] for x in pos]
                   
                    # Draw name
                    self.draw_name(k, list(casAmbiOrph[casAmbiOrph['Operon'] == i]['Prediction'])[0], i, min(starts), max(ends))
                    
                    # Adjust positions
                    startPos =  list(casAmbiOrph[casAmbiOrph['Operon'] == i]['Start'])[0]
                    starts = [self.expand*self.plotexpand + 1 + x - startPos for x in starts]
                    ends = [self.expand*self.plotexpand + 1 + x - startPos for x in ends]

                    cas_list = list(zip(starts, ends, strands, casName, list(((255,0,0),)*len(casName))))
                    
                    # Expand
                    expand_list = self.expandCas(contig, pos, startPos)
                    cas_list = cas_list + expand_list

                    # Draw
                    cas_list = sorted(cas_list, key=lambda x: x[0])
                    self.draw_system(cas_list, [], k)

            # Draw Orphan CRISPR
            if len(self.orphan_crispr) > 0:
                for i in list(self.orphan_crispr['CRISPR']):
                    k += 1
                    logging.debug('Plotting '+i)

                    # Get data
                    contig = list(self.orphan_crispr[self.orphan_crispr['CRISPR'] == i]['Contig'])[0]
                    pred = list(self.orphan_crispr[self.orphan_crispr['CRISPR'] == i]['Prediction'])[0]
                    start = list(self.orphan_crispr[self.orphan_crispr['CRISPR'] == i]['Start'])[0]
                    end = list(self.orphan_crispr[self.orphan_crispr['CRISPR'] == i]['End'])[0]

                    # Expand
                    if self.expand > 0:
                        
                        after_df = self.genes[self.genes['Start'] > end] 
                        before_df = self.genes[self.genes['End'] < start]

                        if len(after_df) > 0:
                            after = after_df.iloc[0,:]['Pos']
                        else:
                            after = 0
                        if len(before_df) > 0:
                            before = before_df.iloc[-1,:]['Pos']
                        else:
                            before = 0

                        # Draw
                        expand_list = self.expandCas(contig, [before, after], start, True)
                        expand_list = sorted(expand_list, key=lambda x: x[0])
                        self.draw_system(expand_list, [], k)
                            

                    # Draw
                    self.draw_array(self.expand*self.plotexpand + 1, self.expand*self.plotexpand + 1 + end - start, pred, k, 1)
                    self.draw_name(k, pred, i, start, end)
                            
            self.im.save(self.out+'plot.png')
                    
