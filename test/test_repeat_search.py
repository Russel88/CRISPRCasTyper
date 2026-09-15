"""Unit tests for the repeat search backends (myers-batch prefilter + edlib)."""
import pytest

from cctyper.editdist import search_repeats_edlib, search_repeats_myers

BACKENDS = [search_repeats_edlib, search_repeats_myers]

# 32 bp repeat -> k = ceil(0.15*32) = 5
REP = 'GTTTCAGACGAACCCTTGTGGGATTGAAGCTC'
RC = 'GAGCTTCAATCCCACAAGGGTTCGTCTGAAAC'
PAD = 'ATATATATCGCGCGCGATATATATCGCGCGCG'  # repeat-free padding


def hits(search, flank, repeats=None):
    rows = search(repeats or [('r1', REP, RC)], {'f1': flank}, threads=2)
    return sorted((r[8], r[9], r[4]) for r in rows)  # (start, end, n_edits)


def mutate(seq, *pos, base='A'):
    s = list(seq)
    for p in pos:
        s[p] = 'T' if s[p] == base else base
    return ''.join(s)


@pytest.mark.parametrize('search', BACKENDS)
def test_exact_and_position(search):
    h = hits(search, PAD + REP + PAD)
    assert (len(PAD)+1, len(PAD)+len(REP), 0) in h


@pytest.mark.parametrize('search', BACKENDS)
def test_mismatches(search):
    assert hits(search, PAD + mutate(REP, 5) + PAD)[0][2] == 1
    assert hits(search, PAD + mutate(REP, 5, 10, 15, 20, 25) + PAD)[0][2] == 5  # k boundary
    assert hits(search, PAD + mutate(REP, 3, 7, 11, 15, 19, 23) + PAD) == []    # 6 > k


@pytest.mark.parametrize('search', BACKENDS)
def test_indels(search):
    assert hits(search, PAD + REP[:16] + REP[18:] + PAD)[0][2] == 2   # 2bp deletion
    assert hits(search, PAD + REP[:16] + 'GGG' + REP[16:] + PAD)[0][2] == 3  # 3bp insertion


@pytest.mark.parametrize('search', BACKENDS)
def test_reverse_complement(search):
    h = hits(search, PAD + RC + PAD)
    assert (len(PAD)+1, len(PAD)+len(REP), 0) in h


@pytest.mark.parametrize('search', BACKENDS)
def test_terminal_truncation(search):
    # last 4 bases missing at flank end: 4 edits <= k, found; 8 missing > k, not
    assert hits(search, PAD + REP[:-4]) != []
    assert hits(search, PAD + REP[:-8]) == []


@pytest.mark.parametrize('search', BACKENDS)
def test_multiple_occurrences(search):
    # 3 copies with varying divergence: all found, not just the best
    flank = PAD + REP + PAD + mutate(REP, 5) + PAD + mutate(REP, 5, 10, 15) + PAD
    assert [h[2] for h in hits(search, flank)] == [0, 1, 3]


@pytest.mark.parametrize('search', BACKENDS)
def test_ambiguous_bases(search):
    # N counts as mismatch, does not crash or match spuriously
    assert hits(search, PAD + mutate(REP, 5, base='N') + PAD)[0][2] == 1
    assert hits(search, 'N' * 200) == []


def test_backends_agree_across_windows():
    # hit far into the text (crosses myers window boundaries at 4096)
    flank = PAD * 200 + REP + PAD * 10
    assert hits(search_repeats_edlib, flank) == hits(search_repeats_myers, flank)
