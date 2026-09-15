"""SVG/PNG map plotting."""
import os

import pytest

from conftest import CAIRO_MAX, GENOMES, genome, png_size, run_cctyper


def check_plot(out):
    assert os.path.getsize(os.path.join(out, 'plot.svg')) > 0
    w, h = png_size(os.path.join(out, 'plot.png'))
    assert 0 < w <= CAIRO_MAX and 0 < h <= CAIRO_MAX, (w, h)


@pytest.mark.parametrize('acc', GENOMES)
def test_plot(acc):
    check_plot(run_cctyper(genome(acc)))


def test_plot_expand():
    # large canvas used to exceed cairo's 32767 px surface limit
    check_plot(run_cctyper(genome('NC_017459.1'), '--expand', '20000'))


def test_no_plot():
    out = run_cctyper(genome('NC_017459.1'), '--no_plot')
    assert not os.path.exists(os.path.join(out, 'plot.png'))
