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
        
        self.font = ImageFont.truetype(os.path.join(self.db, 'arial.ttf'), 3*self.scale)
        self.fontB = ImageFont.truetype(os.path.join(self.db, 'arial.ttf'), 5*self.scale)
        self.fontS = ImageFont.truetype(os.path.join(self.db, 'arial.ttf'), 2*self.scale)

    def draw_gene(self, start, end, strand, name, n, z):
        if strand > 0:
            self.draw.polygon(((self.scale/50*start, n*20*self.scale),
                          (self.scale/50*end-5*self.scale, n*20*self.scale),
                          (self.scale/50*end, 2.5*self.scale+n*20*self.scale),
                          (self.scale/50*end-5*self.scale, 5*self.scale+n*20*self.scale),
                          (self.scale/50*start, 5*self.scale+n*20*self.scale)),
                     fill=(255, 0, 0), outline=(255, 255, 255))
        else:
            self.draw.polygon(((self.scale/50*start+5*self.scale, n*20*self.scale),
                          (self.scale/50*end, n*20*self.scale),
                          (self.scale/50*end, 5*self.scale+n*20*self.scale),
                          (self.scale/50*start+5*self.scale, 5*self.scale+n*20*self.scale),
                          (self.scale/50*start, 2.5*self.scale+n*20*self.scale)),
                     fill=(255, 0, 0), outline=(255, 255, 255))
        name = re.sub('_[0-9]*_.*', '', name)
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
        self.draw.text((self.scale/50*start+self.scale/10, n*20*self.scale-3*self.scale), 'CRISPR({})'.format(subtype), (0,0,0), font=self.font)

    def draw_name(self, n, pred, contig, start, end):
        self.draw.text((self.scale/10, n*20*self.scale-10*self.scale), '{}: {}({}-{})'.format(pred, contig, start, end),
                  (0,0,0), font=self.fontB)

    def draw_system(self, cas, crispr, n):
        z = 0
        if len(cas) > 0:
            for i in cas:
                z += 1
                self.draw_gene(i[0], i[1], i[2], i[3], n, z)
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

    def plot(self):

        total = 0

        # Combine orphan and ambiguous cas operons
        if self.any_operon:
            cas_ambi = self.preddf[self.preddf['Prediction'] == 'Ambiguous']
            try:
                casAmbiOrph = pd.concat([self.orphan_cas, cas_ambi])
            except:
                casAmbiOrph = cas_ambi
            total += len(casAmbiOrph)

        try:
            total += len(self.orphan_crispr)
        except:
            pass

        try:
            total += len(self.crispr_cas)
        except:
            pass

        if (not self.noplot) and total > 0:
            
            logging.info('Plotting map of CRISPR-Cas loci')

            width = self.get_longest(self.orphan_crispr, casAmbiOrph, self.crispr_cas, self.crisprsall) 

            self.im = Image.new('RGB', (int(round(self.scale/50*width+self.scale*2)), int(round((total+1)*20*self.scale))), (255, 255, 255))
            self.draw = ImageDraw.Draw(self.im)

            if self.grid:
                # Draw grid
                y_start = 8*self.scale
                y_end = self.im.height
                step_size = int(round(1000*self.scale/50))

                for x in range(step_size, self.im.width, step_size):
                    line = ((x, y_start), (x, y_end))
                    self.draw.line(line, fill=(0,0,0), width=int(self.scale/20))
                    self.draw.text((x-self.scale*4, self.scale*5), str(int(x/(self.scale/50))), (0,0,0), font=self.fontS)

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
                    crisprs = list(self.crispr_cas[self.crispr_cas['Operon'] == i]['CRISPRs'])
                    startsCris = [list(self.crisprsall[self.crisprsall['CRISPR'] == x[0]]['Start'])[0] for x in crisprs]
                    endsCris = [list(self.crisprsall[self.crisprsall['CRISPR'] == x[0]]['End'])[0] for x in crisprs]
                    nameCris = [list(self.crisprsall[self.crisprsall['CRISPR'] == x[0]]['Prediction'])[0] for x in crisprs]

                    # Find start of loci
                    startPos = min(startsCas + startsCris)

                    # Draw name
                    self.draw_name(k, prediction, i, startPos, max(endsCas + endsCris))
                    
                    # Adjust positions
                    startsCas = [1 + x - startPos for x in startsCas]
                    endsCas = [1 + x - startPos for x in endsCas]
                    startsCris = [1 + x - startPos for x in startsCris]
                    endsCris = [1 + x - startPos for x in endsCris]
                    
                    # Draw
                    cas_list = list(zip(startsCas, endsCas, strands, nameCas))
                    cas_list = sorted(cas_list, key=lambda x: x[0])
                    self.draw_system(cas_list, list(zip(startsCris, endsCris, nameCris)), k)

            # Draw Orphan Cas
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
                    starts = [1 + x - list(casAmbiOrph[casAmbiOrph['Operon'] == i]['Start'])[0] for x in starts]
                    ends = [1 + x - list(casAmbiOrph[casAmbiOrph['Operon'] == i]['Start'])[0] for x in ends]

                    # Draw
                    cas_list = list(zip(starts, ends, strands, casName))
                    cas_list = sorted(cas_list, key=lambda x: x[0])
                    
                    self.draw_system(cas_list, [], k)
                    
            # Draw Orphan CRISPR
            if len(self.orphan_crispr) > 0:
                for i in list(self.orphan_crispr['CRISPR']):
                    k += 1
                    logging.debug('Plotting '+i)

                    # Get data
                    pred = list(self.orphan_crispr[self.orphan_crispr['CRISPR'] == i]['Prediction'])[0]
                    start = list(self.orphan_crispr[self.orphan_crispr['CRISPR'] == i]['Start'])[0]
                    end = list(self.orphan_crispr[self.orphan_crispr['CRISPR'] == i]['End'])[0]

                    # Draw
                    self.draw_array(1, end - start, pred, k, 1)
                    self.draw_name(k, pred, i, start, end)
                            
            self.im.save(self.out+'plot.png')
                    
