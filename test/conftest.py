"""Shared helpers for the cctyper test suite.

Test genomes in test/data/ are precomputed CRISPR-Cas loci (operon + arrays,
+-10 kb) cut from NCBI reference genomes; each fasta header records the
original coordinates. To add a genome: run cctyper on the full genome, extract
the loci regions, and add the accession with its expected subtypes to GENOMES.
"""
import atexit
import csv
import functools
import os
import shutil
import struct
import subprocess
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))

# subset fasta in test/data/ -> expected CRISPR-Cas subtypes
GENOMES = {
    'NC_017459.1': {'I-B', 'I-D'},   # Haloquadratum walsbyi C23
    'NC_000913.3': {'I-E'},          # Escherichia coli K-12 MG1655
    'NC_002737.2': {'I-C', 'II-A'},  # Streptococcus pyogenes M1 GAS
    'NC_002976.3': {'III-A'},        # Staphylococcus epidermidis RP62A
    # Issue #65: array overlapping its operon, repeat so divergent that
    # repeatTyper calls III-D while the cas genes say V-F1
    'DAIF01000015.1': {'V-F1'},      # Ca. Pacearchaeota archaeon UBA73
}

CAIRO_MAX = 32767  # cairo max surface dimension in px

OUTDIRS = []
atexit.register(lambda: [shutil.rmtree(d, ignore_errors=True) for d in OUTDIRS])


def genome(acc):
    return os.path.join(HERE, 'data', acc + '_loci.fasta')


@functools.lru_cache(maxsize=None)  # one cctyper run per unique arg set, shared across tests
def run_cctyper(fasta, *extra_args):
    out = tempfile.mkdtemp(prefix='cctyper_test_')
    shutil.rmtree(out)  # cctyper wants a non-existing dir
    OUTDIRS.append(out)
    res = subprocess.run(['cctyper', fasta, out, *extra_args],
                         capture_output=True, text=True)
    assert res.returncode == 0, res.stderr
    assert 'PNG plot failed' not in res.stderr, res.stderr
    return out


def read_tab(out, name):
    with open(os.path.join(out, name)) as f:
        return list(csv.DictReader(f, delimiter='\t'))


def png_size(path):
    with open(path, 'rb') as f:
        header = f.read(24)
    assert header[:8] == b'\x89PNG\r\n\x1a\n', 'not a PNG'
    return struct.unpack('>II', header[16:24])
