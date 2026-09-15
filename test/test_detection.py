"""CRISPR-Cas detection and output files on each test genome."""
import os
import shutil

import pytest

from conftest import GENOMES, genome, read_tab, run_cctyper


@pytest.mark.parametrize('acc', GENOMES)
def test_subtypes(acc):
    out = run_cctyper(genome(acc))
    preds = {row['Prediction'] for row in read_tab(out, 'CRISPR_Cas.tab')}
    assert preds == GENOMES[acc], preds


@pytest.mark.parametrize('acc', GENOMES)
def test_output_files(acc):
    out = run_cctyper(genome(acc))
    for name in ('cas_operons.tab', 'crisprs_all.tab', 'genes.tab',
                 'hmmer.tab', 'crisprs.gff'):
        assert os.path.getsize(os.path.join(out, name)) > 0, name
    assert os.listdir(os.path.join(out, 'spacers'))


def test_degenerate_array_rescue():
    # I-D array with 4 substitutions per repeat: diced misses it (see --skip_blast),
    # only repeat matching against the known-repeat database finds it
    out = run_cctyper(genome('NC_017459.1_ID_degenerate'))
    assert len(read_tab(out, 'CRISPR_Cas.tab')) == 1
    skip = run_cctyper(genome('NC_017459.1_ID_degenerate'), '--skip_blast')
    assert not os.path.exists(os.path.join(skip, 'CRISPR_Cas.tab'))


@pytest.mark.skipif(not shutil.which('blastn'), reason='blastn not installed')
@pytest.mark.parametrize('acc', GENOMES)
def test_backends_equivalent(acc):
    # the native backends must reproduce what BLAST+ used to give
    ref = run_cctyper(genome(acc), '--repeat_search_backend', 'blast')
    for backend in ('myers', 'edlib'):
        out = run_cctyper(genome(acc), '--repeat_search_backend', backend)
        for name in ('CRISPR_Cas.tab', 'crisprs_near_cas.tab', 'cas_operons.tab'):
            assert read_tab(out, name) == read_tab(ref, name), (acc, backend, name)
