import logging
import math

from concurrent.futures import ThreadPoolExecutor

import edlib
import myers_batch

from Bio import SeqIO
from Bio.Seq import Seq

K_FRACTION = 0.15

WINDOW = 4096
OVERLAP = 64


def _all_hits(query, text, k):
    '''
    All occurrences of query in text with edit distance <= k.
    edlib only reports the best-scoring locations, so mask found hits and repeat
    until nothing is left within the threshold.
    '''
    hits = []
    res = edlib.align(query, text, mode='HW', task='locations', k=k)
    if res['editDistance'] < 0:
        return hits
    masked = bytearray(text, 'ascii')
    while res['editDistance'] >= 0:
        for s, e in res['locations']:
            hits.append((s, e, res['editDistance']))
            masked[s:e+1] = b'N' * (e+1-s)
        res = edlib.align(query, bytes(masked).decode(), mode='HW', task='locations', k=k)
    return hits


def load_repeats(repeatdb):
    repeats = []
    for rec in SeqIO.parse(repeatdb, 'fasta'):
        seq = str(rec.seq).upper()
        repeats.append((rec.id, seq, str(Seq(seq).reverse_complement())))
    return repeats


def _row(rid, fid, s, e, dist, qlen):
    return (rid, fid, round(100*(1-dist/qlen), 1), e-s+1,
            dist, 0, 1, qlen, s+1, e+1, 0, qlen-dist)


def search_repeats_edlib(repeats, flanks, threads):
    '''
    Brute-force scan: edlib alignment of every repeat against every flank.
    '''
    def work(chunk):
        rows = []
        for fid, text in flanks.items():
            for rid, fwd, rev in chunk:
                k = math.ceil(len(fwd) * K_FRACTION)
                for query in (fwd, rev):
                    for s, e, dist in _all_hits(query, text, k):
                        rows.append(_row(rid, fid, s, e, dist, len(fwd)))
        return rows

    n = min(threads, len(repeats))
    with ThreadPoolExecutor(n) as ex:
        results = ex.map(work, [repeats[i::n] for i in range(n)])
    return [r for chunk in results for r in chunk]


def search_repeats_myers(repeats, flanks, threads):
    '''
    Two-stage search: myers-batch SIMD prefilter over flank windows (distance
    only), then edlib on positive windows to recover exact coordinates.
    '''
    logging.debug('myers-batch SIMD backend: %s', myers_batch.simd_backend())

    windows = []
    for fid, text in flanks.items():
        step = WINDOW - OVERLAP
        for off in range(0, max(1, len(text) - OVERLAP), step):
            windows.append((fid, off, text[off:off+WINDOW]))
    wtexts = [w[2].encode() for w in windows]

    def work(chunk):
        rows = set()
        for rid, fwd, rev in chunk:
            k = math.ceil(len(fwd) * K_FRACTION)
            for query in (fwd, rev):
                dists = myers_batch.distances(query.encode(), wtexts)
                for wi, d in enumerate(dists):
                    if d <= k:
                        fid, off, text = windows[wi]
                        for s, e, dist in _all_hits(query, text, k):
                            rows.add(_row(rid, fid, off+s, off+e, dist, len(fwd)))
        return rows

    n = min(threads, len(repeats))
    with ThreadPoolExecutor(n) as ex:
        results = ex.map(work, [repeats[i::n] for i in range(n)])
    return [r for chunk in results for r in chunk]
