import logging
import os
import re
import sys
import io
import pathlib
from functools import partial

import pandas as pd
import pyhmmer
import zstandard
from tqdm import tqdm

class HMMER(object):
    
    def __init__(self, obj):
        self.master = obj
        for key, val in vars(obj).items():
            setattr(self, key, val)

    def main_hmm(self):
        # If redo just get the table
        if self.redo:
            self.read_hmm()
        # Else run HMMER load and write data
        else:
            self.run_hmm()
            self.load_hmm()
            self.write_hmm()
            self.run_custom_hmm()

        # Check if any cas genes
        self.check_hmm()

        # Parse
        self.parse_hmm()

        # Load Custom HMM db
        self.load_custom_hmm()

    # A single search
    def run_hmm(self):
        
        logging.info('Running HMMER (pyhmmer) against Cas profiles')

        os.makedirs(self.out+'hmmer', exist_ok=True)

        alphabet = pyhmmer.easel.Alphabet.amino()
        with pyhmmer.easel.SequenceFile(self.prot_path, digital=True, alphabet=alphabet) as seq_file:
            sequences = list(seq_file)
        self.alphabet = alphabet
        self.sequences = sequences

        seq_lookup = {seq.name: seq for seq in sequences}

        hits_out = []
        hmms = list(self._iter_hmms(alphabet))
        cpus = self.threads if self.threads and self.threads > 0 else 0
        progress = partial(
            tqdm,
            total=len(hmms),
            disable=getattr(self, "simplelog", False),
            unit="hmm",
            desc="HMMER",
            leave=False,
        )

        for hits in progress(pyhmmer.hmmsearch(hmms, sequences, cpus=cpus)):
            hmm_name = hits.query.name
            qlen = hits.query.M
            for hit in hits.included:
                target_name = hit.name
                seq = seq_lookup.get(target_name)
                if seq is None:
                    continue
                start = end = strand = 0  # gene table will fill these later
                tlen = len(seq)
                for domain in hit.domains.included:
                    aln = domain.alignment
                    if aln is not None:
                        t_from = aln.target_from
                        t_to = aln.target_to
                        h_from = aln.hmm_from
                        h_to = aln.hmm_to
                    else:
                        t_from = domain.env_from
                        t_to = domain.env_to
                        h_from = 0
                        h_to = 0

                    cov_seq = (t_to - t_from + 1) / tlen if tlen else 0
                    if qlen:
                        if h_to and h_from:
                            cov_hmm = (h_to - h_from + 1) / qlen
                        else:
                            cov_hmm = (domain.env_to - domain.env_from + 1) / qlen
                    else:
                        cov_hmm = 0
                    if cov_seq <= 0:
                        cov_seq = 1
                    if cov_hmm <= 0:
                        cov_hmm = 1
                    hits_out.append({
                        "Hmm": hmm_name,
                        "ORF": target_name,
                        "tlen": tlen,
                        "qlen": qlen,
                        "Eval": domain.i_evalue,
                        "score": domain.score,
                        "hmm_from": h_from,
                        "hmm_to": h_to,
                        "ali_from": t_from,
                        "ali_to": t_to,
                        "env_from": domain.env_from,
                        "env_to": domain.env_to,
                        "pprop": getattr(domain, "accuracy", None),
                        "start": start,
                        "end": end,
                        "strand": strand,
                        "Cov_seq": cov_seq,
                        "Cov_hmm": cov_hmm,
                    })

        if not hits_out:
            self.hmm_df = pd.DataFrame(columns=['Hmm','ORF','tlen','qlen','Eval','score','start','end','Acc','Pos','Cov_seq','Cov_hmm','strand'])
            return

        hmm_df = pd.DataFrame(hits_out)

        genes = pd.read_csv(self.out+'genes.tab', sep='\t')

        # Drop placeholder coords so gene-derived values don't create duplicates
        hmm_df = hmm_df.drop(columns=['start', 'end', 'strand'], errors='ignore')

        if getattr(self, "gff", False) and getattr(self, "prot", False):
            hmm_df = pd.merge(
                hmm_df,
                genes[["Start", "End", "Strand", "Contig", "Pos", "protein_id"]],
                left_on="ORF",
                right_on="protein_id",
                how="left",
            ).drop("protein_id", axis=1)
            hmm_df.rename(columns={'Contig': 'Acc','Strand': 'strand', 'Start': 'start', 'End': 'end'}, inplace=True)
        else:
            def parse_acc_pos(name: str):
                if '#' in name:
                    core = name.split('#')[0].strip()
                    return re.sub("_[0-9]*$","", core), int(re.sub(".*_","", core))
                m = re.search(r"_(\d+)$", name)
                if m:
                    return re.sub("_[0-9]*$","",name), int(m.group(1))
                return name, 0

            acc_pos = [parse_acc_pos(x) for x in hmm_df['ORF']]
            hmm_df['Acc'] = [ap[0] for ap in acc_pos]
            hmm_df['Pos'] = [ap[1] for ap in acc_pos]

            hmm_df = pd.merge(
                hmm_df,
                genes[["Start", "End", "Strand", "Contig", "Pos"]],
                left_on=["Acc","Pos"],
                right_on=["Contig","Pos"],
                how="left",
            ).drop("Contig", axis=1)

        # Prefer gene-derived coordinates/strand over header defaults (zeros)
        for src, tgt in (('Start', 'start'), ('End', 'end'), ('Strand', 'strand')):
            if src in hmm_df.columns:
                hmm_df[tgt] = pd.to_numeric(hmm_df[src], errors='coerce')

        # Drop merge artefacts and duplicates
        for col in ['Start', 'End', 'Strand', 'start_x', 'end_x', 'strand_x', 'start_y', 'end_y', 'strand_y']:
            if col in hmm_df.columns:
                hmm_df.drop(columns=[col], inplace=True)

        hmm_df = hmm_df.loc[:, ~hmm_df.columns.duplicated()]
        self.hmm_df = hmm_df.drop_duplicates()

    def _iter_hmms(self, alphabet):

        def load_plain(path: pathlib.Path):
            with pyhmmer.plan7.HMMFile(path) as hmm_file:
                for hmm in hmm_file:
                    yield hmm

        if self.hmm_db.endswith('.zst'):
            target_plain = pathlib.Path(self.hmm_db).with_suffix('')
            if target_plain.exists():
                try:
                    yield from load_plain(target_plain)
                    return
                except Exception as exc:
                    logging.warning('Existing decompressed HMM %s unreadable (%s); regenerating from zst', target_plain, exc)

            with open(self.hmm_db, 'rb') as fh:
                compressed = fh.read()
            dctx = zstandard.ZstdDecompressor()
            try:
                data = dctx.decompress(compressed)
            except zstandard.ZstdError as exc:
                logging.error('Failed to decompress %s: %s', self.hmm_db, exc)
                sys.exit()

            target_plain.write_bytes(data)
            yield from load_plain(target_plain)
        else:
            yield from load_plain(pathlib.Path(self.hmm_db))

    def _parse_prodigal_header(self, seq):
        desc = (seq.description or "")
        start = end = strand = 0
        if '#' in desc:
            parts = [p.strip() for p in desc.split('#')]
            # Headers we generate look like:
            #   >contig_idx # start # end # strand               (prodigal)
            #   >pid # contig_pos # start # end # strand         (gff/prot)
            # Always take the last three numeric fields if present
            if len(parts) >= 4:
                try:
                    start = int(parts[-3])
                    end = int(parts[-2])
                    strand = int(parts[-1])
                except ValueError:
                    start = end = strand = 0
        return start, end, strand

    # Load data
    def load_hmm(self):
    
        logging.debug('Loading HMMER output')
        
        def merge_hits(df_sub):
            tlen = df_sub['tlen'].iloc[0]
            qlen = df_sub['qlen'].iloc[0]

            start_series = df_sub['start']
            end_series = df_sub['end']
            strand_series = df_sub['strand'] if 'strand' in df_sub.columns else None
            if isinstance(start_series, pd.DataFrame):
                start_series = start_series.iloc[:, 0]
            if isinstance(end_series, pd.DataFrame):
                end_series = end_series.iloc[:, 0]
            if isinstance(strand_series, pd.DataFrame):
                strand_series = strand_series.iloc[:, 0]
            start_series = pd.to_numeric(start_series, errors='coerce')
            end_series = pd.to_numeric(end_series, errors='coerce')
            if strand_series is not None:
                strand_series = pd.to_numeric(strand_series, errors='coerce')

            seq_span = set()
            for i, j in zip(df_sub['ali_from'], df_sub['ali_to']):
                seq_span.update(range(int(i), int(j) + 1))
            hmm_span = set()
            for i, j in zip(df_sub['hmm_from'], df_sub['hmm_to']):
                if i and j:
                    hmm_span.update(range(int(i), int(j) + 1))
            cov_seq = len(seq_span) / tlen if tlen else 0
            cov_hmm = len(hmm_span) / qlen if qlen else 0

            best = df_sub.sort_values('score', ascending=False).iloc[0].copy()
            valid_start = start_series[start_series > 0]
            valid_end = end_series[end_series > 0]
            best['start'] = valid_start.min() if not valid_start.empty else start_series.min()
            best['end'] = valid_end.max() if not valid_end.empty else end_series.max()
            if strand_series is not None:
                strands = strand_series[strand_series != 0].dropna()
                if not strands.empty:
                    best['strand'] = strands.iloc[0]
            best['Cov_seq'] = cov_seq
            best['Cov_hmm'] = cov_hmm
            return best

        def merge_hits_with_keys(df_sub):
            merged = merge_hits(df_sub)
            if 'Hmm' not in merged or 'ORF' not in merged:
                if hasattr(df_sub, "name"):
                    hmm_key, orf_key = df_sub.name
                    merged['Hmm'] = hmm_key
                    merged['ORF'] = orf_key
            return merged

        self.hmm_df = (
            self.hmm_df
            .groupby(['Hmm','ORF'], group_keys=False)
            .apply(merge_hits_with_keys, include_groups=False)
            .reset_index(drop=True)
        )

        cols = ['Hmm','ORF','tlen','qlen','Eval','score','start','end','Acc','Pos','Cov_seq','Cov_hmm','strand']
        self.hmm_df = self.hmm_df[cols]

    # Write to file
    def write_hmm(self):
        self.hmm_df.to_csv(self.out+'hmmer.tab', sep='\t', index=False)

    # Read from file
    def read_hmm(self):
        try:        
            self.hmm_df = pd.read_csv(self.out+'hmmer.tab', sep='\t')
        except:
            logging.error('No matches to Cas HMMs')
            sys.exit()

    # Check if any cas genes
    def check_hmm(self):
        if len(self.hmm_df) == 0:
            logging.info('No Cas proteins found.')
        else:
            self.any_cas = True

    # Parse
    def parse_hmm(self):

        if self.any_cas:
        
            logging.debug('Parsing HMMER output')

            # Pick best hit
            self.hmm_df.sort_values('score', ascending=False, inplace=True)
            self.hmm_df.drop_duplicates('ORF', inplace=True)

    def run_custom_hmm(self):
        
        if self.customhmm != '':
            logging.info('Running HMMER (pyhmmer) against custom HMM profiles')

            if not hasattr(self, 'sequences'):
                alphabet = pyhmmer.easel.Alphabet.amino()
                with pyhmmer.easel.SequenceFile(self.prot_path, digital=True, alphabet=alphabet) as seq_file:
                    self.sequences = list(seq_file)
                self.alphabet = alphabet

            hits_out = []
            cpus = self.threads if self.threads and self.threads > 0 else 0

            with pyhmmer.plan7.HMMFile(self.customhmm) as hmm_file:
                hmms = list(hmm_file)

            for hits in pyhmmer.hmmsearch(hmms, self.sequences, cpus=cpus):
                hmm_name = (hits.query.name or "")
                hmm_acc = (hits.query.accession or "")
                for hit in hits.included:
                    target_name = hit.name
                    for domain in hit.domains.included:
                        aln = domain.alignment
                        if aln is not None:
                            t_from = aln.target_from
                            t_to = aln.target_to
                            h_from = aln.hmm_from
                            h_to = aln.hmm_to
                        else:
                            t_from = domain.env_from
                            t_to = domain.env_to
                            h_from = 0
                            h_to = 0
                        hits_out.append({
                            "Target": target_name,
                            "Query": hmm_name,
                            "Acc": hmm_acc,
                            "E-value": domain.i_evalue,
                            "Score": domain.score,
                            "ali_from": t_from,
                            "ali_to": t_to,
                            "hmm_from": h_from,
                            "hmm_to": h_to,
                        })

            if hits_out:
                self.custom_hmm_df = pd.DataFrame(hits_out)
            else:
                self.custom_hmm_df = pd.DataFrame(columns=['Target','Query','Acc','E-value','Score'])
            
    def load_custom_hmm(self):

        if self.customhmm != '':

            # Remove low E-value hits
            self.custom_hmm_df = self.custom_hmm_df[self.custom_hmm_df['E-value'] < self.oev]
            
            # Pick best hit
            self.custom_hmm_df.sort_values('Score', ascending=False, inplace=True)
            self.custom_hmm_df.drop_duplicates('Target', inplace=True)

            # New columns
            self.custom_hmm_df['Contig'] = [re.sub("_[0-9]*$","",x) for x in self.custom_hmm_df['Target']]
            self.custom_hmm_df['Pos'] = [int(re.sub(".*_","",x)) for x in self.custom_hmm_df['Target']]
