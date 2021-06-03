import os
import logging
import re

import pandas as pd

import drawSvg as draw

class Map(object):
    
    def __init__(self, obj):
        self.master = obj
        for key, val in vars(obj).items():
            setattr(self, key, val)
        
    def draw_gene(self, start, end, strand, name, n, z, put):

        if isinstance(name, str):
            font_size = 26
            name_full = name
            name = re.sub('_[0-9]*_.*', '', name)
            if name_full in self.cas_hmms:
                if 'Cas6' in name:
                    col = '#dd1424'
                elif 'Cas3-Cas2' in name:
                    col = '#009b77'
                elif any([True if x in name else False for x in self.interf_genes]):
                    col = '#efc050'
                elif any([True if x in name+'_' else False for x in self.adapt_genes]):
                    col = '#5b5ea6'
                else:
                    col = '#ff80e6'
            else:
                col = '#d3c9c1'

        else:
            font_size = 20
            name = str(name)
            col = '#d3c9c1'

        opa = 0.8
        if put:
            name = '('+name+')'
            opa = 0.3

        if strand > 0:
            self.im.append(draw.Lines(
                10+self.scale/50*start, self.imheight-(n*20*self.scale),
                10+self.scale/50*end-3*self.scale, self.imheight-(n*20*self.scale),
                10+self.scale/50*end, self.imheight-(2.5*self.scale+n*20*self.scale),
                10+self.scale/50*end-3*self.scale, self.imheight-(5*self.scale+n*20*self.scale),
                10+self.scale/50*start, self.imheight-(5*self.scale+n*20*self.scale),
                10+self.scale/50*start, self.imheight-(n*20*self.scale),
                     fill=col, close=True, fill_opacity=opa, stroke=col, stroke_width=4))
        else:
            self.im.append(draw.Lines(
                10+self.scale/50*start+3*self.scale, self.imheight-(n*20*self.scale),
                10+self.scale/50*end, self.imheight-(n*20*self.scale),
                10+self.scale/50*end, self.imheight-(5*self.scale+n*20*self.scale),
                10+self.scale/50*start+3*self.scale, self.imheight-(5*self.scale+n*20*self.scale),
                10+self.scale/50*start, self.imheight-(2.5*self.scale+n*20*self.scale),
                10+self.scale/50*start+3*self.scale, self.imheight-(n*20*self.scale),
                     fill=col, close=True, fill_opacity=opa, stroke=col, stroke_width=4))
        
        if z % 2 == 1:
            self.im.append(draw.Text(name, font_size, self.scale/50*start+10, self.imheight-(n*20*self.scale-1*self.scale), fill='black'))
        else:
            self.im.append(draw.Text(name, font_size, self.scale/50*start+10, self.imheight-(n*20*self.scale+8*self.scale), fill='black'))

    def draw_array(self, start, end, subtype, n, z, n_reps):
        self.im.append(draw.Lines(
                    10+self.scale/50*start, self.imheight-(n*20*self.scale), 
                    10+self.scale/50*end, self.imheight-(n*20*self.scale), 
                    10+self.scale/50*end, self.imheight-(5*self.scale+n*20*self.scale), 
                    10+self.scale/50*start, self.imheight-(5*self.scale+n*20*self.scale),
                    10+self.scale/50*start, self.imheight-(n*20*self.scale),
                     close=True, fill='white', stroke='white', stroke_width=4))
        
        obj_width = ((end-start)*self.scale/50) / ((n_reps*2)-1)
        
        y_start = self.imheight-(n*20*self.scale)
        y_end = self.imheight-(5*self.scale+n*20*self.scale)

        x_start = 10+self.scale/50*start
        x_end = x_start + obj_width

        obj_type = 1
        for cris_obj in range(((n_reps*2)-1)):
            obj_type += 1
            if obj_type % 2 == 1:
                self.im.append(draw.Lines(
                            x_start, y_start, x_end, y_start, x_end, y_end, x_start, y_end, x_start, y_start,
                            close=True, fill='#f0f0f0'))
            else:
                self.im.append(draw.Lines(
                            x_start, y_start, x_end, y_start, x_end, y_end, x_start, y_end, x_start, y_start,
                            close=True, fill='black'))
            x_start += obj_width
            x_end += obj_width

        if subtype == 'Unknown':
            self.im.append(draw.Text('CRISPR', 26, 10+self.scale/50*start, self.imheight-(n*20*self.scale-1*self.scale), fill='black'))
        else:
            self.im.append(draw.Text('CRISPR: '+subtype, 26, 10+self.scale/50*start, self.imheight-(n*20*self.scale-1*self.scale), fill='black'))

    def draw_name(self, n, pred, contig, start, end):
        self.im.append(draw.Text('{}: {} ({}-{})'.format(pred, contig, start, end), 38, 15+self.scale/10, self.imheight-(n*20*self.scale-6*self.scale), fill='black'))

    def draw_system(self, cas, crispr, n):
        z = 0
        if len(cas) > 0:
            for i in cas:
                z += 1
                self.draw_gene(i[0], i[1], i[2], i[3], n, z, i[4])
        if len(crispr) > 0:
            for i in crispr:
                z += 1
                self.draw_array(i[0], i[1], i[2], n, z, i[3])

    def criscas_len(self, cc, cca):
        
        lengths = []
        self.criscaspos = {}
        for index, row in cc.iterrows():
            span_ends = False
            start = row['Operon_Pos'][0] 
            end = row['Operon_Pos'][1]
            # If Cas operon span ends of sequences
            if start > end:
                span_ends = True
            # If CRISPR-Cas does not span ends
            if not (row['Operon'] in self.cc_circ_start.keys() or row['Operon'] in self.cc_circ_end.keys()):
                start = min(start, min(cca[cca['CRISPR'].isin(row['CRISPRs'])]['Start']))
                end = max(end, max(cca[cca['CRISPR'].isin(row['CRISPRs'])]['End']))
            # If CRISPR-Cas span ends and array is in start of sequence
            if row['Operon'] in self.cc_circ_start.keys():
                ccs = cca[cca['CRISPR'].isin(self.cc_circ_start[row['Operon']])]
                end = max(ccs['End'])
                span_ends = True
            # If CRISPR-Cas span ends and array is in end of sequence
            if row['Operon'] in self.cc_circ_end.keys():
                ccs = cca[cca['CRISPR'].isin(self.cc_circ_end[row['Operon']])]
                start = min(ccs['Start'])
                span_ends = True
            
            seq_size = self.len_dict[row['Contig']]
            
            if span_ends:
                lengths.append(end + seq_size - start)
            else:
                lengths.append(end - start)
            
            self.criscaspos[row['Operon']] = [start, end, seq_size, span_ends]

        return lengths
            
    def get_longest(self, crisO, casO, criscasO, crisA):
        crisO_M, casO_M, cc_M = 0, 0, 0
        if len(crisO) > 0:
            crisO_M = max(crisO['End']-crisO['Start'])
        if len(casO) > 0:
            casO_lin = casO[casO['Start'] < casO['End']]
            casO_circ = casO[casO['Start'] > casO['End']]
            casO_lin_M, casO_circ_M = 0, 0
            if len(casO_lin) > 0:
                casO_lin_M = max(casO_lin['End']-casO_lin['Start'])
            if len(casO_circ) > 0:
                casO_circ['size'] = casO_circ.apply(lambda x: self.len_dict[x['Contig']], axis=1)
                casO_circ_M = max(casO_circ['Start']+casO_circ['size']-casO_circ['Start'])
            casO_M = max(casO_lin_M, casO_circ_M)
        if len(criscasO) > 0:
            cc_M = max(self.criscas_len(criscasO, crisA))
        return(max(crisO_M, casO_M, cc_M))

    def expandCas(self, contig, pos, startPos, endPos, seq_size, span_ends, array=False):
        all_cas = list(self.genes[self.genes['Contig'] == contig]['Pos'])
        
        if array:
            missing_cas = all_cas
        else:
            missing_cas = [x for x in all_cas if x not in pos]
        
        add_these = self.genes[(self.genes['Contig'] == contig) & (self.genes['Pos'].isin(missing_cas))]
        
        if span_ends:
            add_these = add_these[(add_these['Start'] > startPos-self.expand) | (add_these['End'] < endPos+self.expand)]
            which_end = [x>(endPos+self.expand) for x in list(add_these['Start'])]
            add_starts = [self.expand + 1 + x[0] - startPos if x[1] else self.expand + 1 + x[0] + seq_size - startPos for x in zip(list(add_these['Start']),which_end)]
            add_ends = [self.expand + 1 + x[0] - startPos if x[1] else self.expand + 1 + x[0] + seq_size - startPos for x in zip(list(add_these['End']),which_end)]
        else:
            add_these = add_these[(add_these['Start'] > startPos-self.expand) & (add_these['End'] < endPos+self.expand)]
            add_starts = [self.expand + 1 + x - startPos for x in list(add_these['Start'])]
            add_ends = [self.expand + 1 + x - startPos for x in list(add_these['End'])]

        names = list(add_these['Pos'])
        puts = list((False,)*len(add_starts))

        if self.customhmm != '' and len(self.custom_hmm_df)>0: 
            # Add custom
            hmm_contig = self.custom_hmm_df[self.custom_hmm_df['Contig'] == contig]
            add_custom = [x in list(hmm_contig['Pos']) for x in list(add_these['Pos'])]
            custom_names = [list(hmm_contig[hmm_contig['Pos'] == x]['Query']) for x in list(add_these['Pos'])]
            custom_names = [x[0] if len(x)>0 else x for x in custom_names]
            names = [x[0] if not x[1] else x[2] for x in zip(names, add_custom, custom_names)]
        else:
            add_custom = puts

        # Add putative
        try:
            hmm_contig = self.hmm_df_raw[self.hmm_df_raw['Acc'] == contig]
            add_putative = [x in list(hmm_contig['Pos']) for x in list(add_these['Pos'])]
            add_putative = [x[0] if not x[1] else False for x in zip(add_putative, add_custom)]
            put_names = [list(hmm_contig[hmm_contig['Pos'] == x]['Hmm']) for x in list(add_these['Pos'])]
            put_names = [x[0] if len(x)>0 else x for x in put_names]
            puts = [x[0] if not x[1] else True for x in zip(puts, add_putative)]
            names = [x[0] if not x[1] else x[2] for x in zip(names, add_putative, put_names)]
        except:
            pass

        expand_list = list(zip(add_starts,
                          add_ends,
                          list(add_these['Strand']),
                          names,
                          puts))
        
        return expand_list    

    def expandCris(self, contig, crisprs, startPos, endPos, seq_size, span_ends):
      
        add_crisp = self.crisprsall[self.crisprsall['Contig'] == contig]
        add_crisp = add_crisp[~add_crisp['CRISPR'].isin(crisprs)]
        
        if span_ends:
            add_crisp = add_crisp[(add_crisp['Start'] > startPos-self.expand) | (add_crisp['End'] < endPos+self.expand)]
        else:
            add_crisp = add_crisp[(add_crisp['Start'] > startPos-self.expand) & (add_crisp['End'] < endPos+self.expand)]

        if len(add_crisp) > 0:
            crisp_lst = list(add_crisp['CRISPR'])
            startsCris = [list(add_crisp[add_crisp['CRISPR'] == x]['Start'])[0] for x in crisp_lst]
            endsCris = [list(add_crisp[add_crisp['CRISPR'] == x]['End'])[0] for x in crisp_lst]
            nameCris = [list(add_crisp[add_crisp['CRISPR'] == x]['Prediction'])[0] for x in crisp_lst]
            repsCris = [list(add_crisp[add_crisp['CRISPR'] == x]['N_repeats'])[0] for x in crisp_lst]
            
            if span_ends:
                which_end = [x>(endPos+self.expand) for x in startsCris]
                startsCris = [self.expand + 1 + x[0] - startPos if x[1] else self.expand + 1 + x[0] + seq_size - startPos for x in zip(startsCris,which_end)]
                endsCris = [self.expand + 1 + x[0] - startPos if x[1] else self.expand + 1 + x[0] + seq_size - startPos for x in zip(endsCris,which_end)]
            else:
                startsCris = [self.expand + 1 + x - startPos for x in startsCris]
                endsCris = [self.expand + 1 + x - startPos for x in endsCris]

            return list(zip(startsCris, endsCris, nameCris, repsCris))
        else:
            return []

    def plot(self):

        total = 0

        # Combine orphan and ambiguous cas operons
        if self.any_operon:
            cas_ambi = self.preddf[self.preddf['Prediction'] == 'Ambiguous']
            try:
                casAmbiOrph = pd.concat([self.orphan_cas, cas_ambi])
                # Remove ambiguous, which are in CRISPR-Cas
                casAmbiOrph = casAmbiOrph[~casAmbiOrph['Operon'].isin(self.crispr_cas['Operon'])]
            except:
                cas_good = self.preddf[~self.preddf['Prediction'].isin(['False','Ambiguous','Partial'])]
                casAmbiOrph = pd.concat([cas_good, cas_ambi])
            total += len(casAmbiOrph)
        else:
            casAmbiOrph = []
        
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
            
            self.interf_genes = set([subitem for subsublist in [item for sublist in self.compl_interf.values() for item in sublist] for subitem in subsublist])
            self.adapt_genes = set([subitem for subsublist in [item for sublist in self.compl_adapt.values() for item in sublist] for subitem in subsublist])
        
            width = self.get_longest(self.orphan_crispr, casAmbiOrph, self.crispr_cas, self.crisprsall) 

            self.genes = pd.read_csv(self.out+'genes.tab', sep='\t') 
           
            # Set scale
            self.scale = 10

            # Make empty canvas
            width = width + (self.expand * 2)
            self.imheight = int(round((total+1)*20*self.scale))

            self.im = draw.Drawing(int(round(self.scale/50*width+self.scale*10)), self.imheight, displayInline=False, origin=(0,0))

            if not self.nogrid:
                # Draw grid
                y_start = 20
                y_end = self.imheight-8*self.scale
                step_size = int(round(1000*self.scale/50))

                for x in range(step_size, self.im.width, step_size):
                    self.im.append(draw.Lines(x+10, y_start, x+10, y_end, stroke='grey', stroke_width=1, fill='none'))
                    self.im.append(draw.Text(str(int(x/(self.scale/50))), 22, x-40, self.imheight-7*self.scale, fill='grey'))

            # Init count of loci
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
                    repsCris = [list(self.crisprsall[self.crisprsall['CRISPR'] == x]['N_repeats'])[0] for x in crisprs]

                    # Find start of loci
                    startPos = self.criscaspos[i][0]
                    endPos = self.criscaspos[i][1]

                    # Draw name
                    self.draw_name(k, prediction, i, startPos, endPos)
                    
                    # Adjust positions
                    seq_size = self.len_dict[contig]
                    
                    if self.criscaspos[i][3]:
                        startsCas = [self.expand + 1 + x - startPos if x>=startPos else self.expand + 1 + x + seq_size - startPos for x in startsCas]
                        endsCas = [self.expand + 1 + x - startPos if x>startPos else self.expand + 1 + x + seq_size - startPos for x in endsCas]
                        startsCris = [self.expand + 1 + x - startPos if x>=startPos else self.expand + 1 + x + seq_size - startPos for x in startsCris]
                        endsCris = [self.expand + 1 + x - startPos if x>startPos else self.expand + 1 + x + seq_size - startPos for x in endsCris]
                        self.im.append(draw.Lines(
                               10+(self.expand+1+seq_size-startPos)*self.scale/50, self.imheight-(k*20*self.scale-5*self.scale),
                               10+(self.expand+1+seq_size-startPos)*self.scale/50, self.imheight-(k*20*self.scale+10*self.scale),
                               stroke='black', stroke_width=5, fill='none'))
                        
                    else:
                        startsCas = [self.expand + 1 + x - startPos for x in startsCas]
                        endsCas = [self.expand + 1 + x - startPos for x in endsCas]
                        startsCris = [self.expand + 1 + x - startPos for x in startsCris]
                        endsCris = [self.expand + 1 + x - startPos for x in endsCris]
                        
                    # Draw
                    cas_list = list(zip(startsCas, endsCas, strands, nameCas, list((False,)*len(nameCas))))

                    # Expand
                    expand_list = self.expandCas(contig, posCas, startPos, endPos, seq_size, self.criscaspos[i][3])
                    cas_list = cas_list + expand_list
                    expand_cris = self.expandCris(contig, crisprs, startPos, endPos, seq_size, self.criscaspos[i][3])
                    cris_list = list(zip(startsCris, endsCris, nameCris, repsCris)) + expand_cris

                    cas_list = sorted(cas_list, key=lambda x: x[0])
                    self.draw_system(cas_list, cris_list, k)
            
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
                    endPos =  list(casAmbiOrph[casAmbiOrph['Operon'] == i]['End'])[0]
                    seq_size = self.len_dict[contig]
                    
                    if startPos < endPos:
                        span_ends = False
                        starts = [self.expand + 1 + x - startPos for x in starts]
                        ends = [self.expand + 1 + x - startPos for x in ends]
                    else:
                        span_ends = True
                        starts = [self.expand + 1 + x - startPos if x>=startPos else self.expand + 1 + x + seq_size - startPos for x in starts]
                        ends = [self.expand + 1 + x - startPos if x>startPos else self.expand + 1 + x + seq_size - startPos for x in ends]
                        self.im.append(draw.Lines(
                               10+(self.expand+1+seq_size-startPos)*self.scale/50, self.imheight-(self.scale*k*20-5*self.scale),
                               10+(self.expand+1+seq_size-startPos)*self.scale/50, self.imheight-(self.scale*k*20+10*self.scale),
                               stroke='black', stroke_width=5, fill='none'))

                    cas_list = list(zip(starts, ends, strands, casName, list((False,)*len(casName))))
                    
                    # Expand
                    expand_list = self.expandCas(contig, pos, startPos, endPos, seq_size, span_ends)
                    cas_list = cas_list + expand_list
                    
                    if self.any_crispr:
                        expand_cris = self.expandCris(contig, [], startPos, endPos, seq_size, span_ends)
                    else:
                        expand_cris = []

                    # Draw
                    cas_list = sorted(cas_list, key=lambda x: x[0])
                    self.draw_system(cas_list, expand_cris, k)

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
                    reps = list(self.orphan_crispr[self.orphan_crispr['CRISPR'] == i]['N_repeats'])[0]

                    # Expand
                    if self.expand > 0:
                        
                        # Draw
                        expand_list = self.expandCas(contig, [0], start, end, 0, False, True)
                        expand_list = sorted(expand_list, key=lambda x: x[0])
                        expand_cris = self.expandCris(contig, [i], start, end, 0, False)
                        self.draw_system(expand_list, expand_cris, k)
                            
                    # Draw
                    self.draw_array(self.expand + 1, self.expand + 1 + end - start, pred, k, 1, reps)
                    self.draw_name(k, pred, i, start, end)
                    
            self.im.saveSvg(self.out+'plot.svg')
            try: 
                self.im.setPixelScale(int(round(self.im.width/(250*self.scale))))
                self.im.savePng(self.out+'plot.png')
            except:
                logging.warning('PNG plot failed. Trying lower resolution')
                try:
                    self.im.setPixelScale(3)
                    self.im.savePng(self.out+'plot.png')
                except:
                    logging.warning('PNG plot failed')
