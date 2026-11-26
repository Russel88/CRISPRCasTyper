"""An input without Cas genes is a valid result, not a crash.

cctyper used to abort such runs twice: load_hmm() groupby'ed the empty hit table
and lost its columns (KeyError: "['Hmm', 'ORF'] not in index"), and typing()
assumed at least one hit (UnboundLocalError: 'specifics'). Both raised although
"no CRISPR-Cas found" is a valid answer, so genomes without a CRISPR-Cas system
could not be analysed.
"""
import os

from conftest import read_tab, run_cctyper

HERE = os.path.dirname(os.path.abspath(__file__))
# 6 kb of random sequence: pyrodigal finds genes, HMMER finds no Cas protein
NO_CAS = os.path.join(HERE, 'data', 'no_cas_6kb.fasta')


def test_no_cas_hits():
    out = run_cctyper(NO_CAS, '--prodigal', 'meta')
    # no hits, no operons, but a complete run
    assert read_tab(out, 'hmmer.tab') == []
    assert not os.path.exists(os.path.join(out, 'CRISPR_Cas.tab'))
    assert not os.path.exists(os.path.join(out, 'cas_operons.tab'))
    # the cleanup step was reached
    assert not os.path.exists(os.path.join(out, 'proteins.faa'))
