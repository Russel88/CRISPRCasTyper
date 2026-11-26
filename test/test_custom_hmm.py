"""Runs with --custom_hmm must finish, cleanup included.

cctyper used to end every run that used --custom_hmm with a FileNotFoundError:
clean() still stat'ed 'hmmer_custom.log', the log of the old, subprocess based
custom search, which the pyhmmer based search does not write. All results,
including the plot, were written before that, but the temporary files were left
behind and the process exited non-zero.
"""
import os

from conftest import genome, run_cctyper

# HMMER3 model built from one protein of the NC_000913.3 test locus (the
# labels cctyper shows for it are not used by the tests)
CUSTOM_HMM = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data', 'custom.hmm')


def test_custom_hmm_run():
    # a plain run: this used to fail in clean(), after all results were written
    out = run_cctyper(genome('NC_000913.3'), '--custom_hmm', CUSTOM_HMM)
    assert os.path.exists(os.path.join(out, 'CRISPR_Cas.tab'))
    # the cleanup step ran, so the temporary files are gone
    assert not os.path.exists(os.path.join(out, 'proteins.faa'))
    assert not os.path.exists(os.path.join(out, 'hmmer'))
    assert not os.path.exists(os.path.join(out, 'hmmer_custom.log'))


def test_custom_hmm_run_without_plot():
    # the crash was independent of the plot, which the custom profiles decorate
    out = run_cctyper(genome('NC_000913.3'), '--custom_hmm', CUSTOM_HMM, '--no_plot')
    assert not os.path.exists(os.path.join(out, 'plot.svg'))
    assert not os.path.exists(os.path.join(out, 'proteins.faa'))
