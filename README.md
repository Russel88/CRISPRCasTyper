# CasPredict

Detect CRISPR-Cas genes, group them into operons, and predict their subtype 

## Installation
### Conda
It is advised to use [miniconda](https://docs.conda.io/en/latest/miniconda.html) or [anaconda](https://www.anaconda.com/) to install.

```sh
conda create -n caspredict -c russel88 caspredict
```

### pip
However, if you have the dependencies (Python >= 3.8, HMMER >= 3.2, Prodigal >= 2.6, grep, sed) in your PATH you can install with pip

```sh
python -m pip install caspredict
```

## Download database
### Conda
Coming soon...

### pip
Coming soon...

## How to run
##### Activate environment
```sh
conda activate caspredict
```
##### Run with a nucleotide fasta as input
```sh
caspredict genome.fa my_output
```
##### Use multiple threads
```sh
caspredict genome.fa my_output -t 20
```
##### Check the different options
```sh
caspredict -h

usage: caspredict [-h] [-t THREADS] [--prodigal {single,meta}] [--aa] [--skip_check] [--keep_prodigal] [--log_lvl {DEBUG,INFO,WARNING,ERROR}] [--redo_typing] [--scores SCORES] [--hmms HMMS]
                  [--dist DIST] [--overall_eval OVERALL_EVAL] [--overall_cov_seq OVERALL_COV_SEQ] [--overall_cov_hmm OVERALL_COV_HMM] [--two_gene_eval TWO_GENE_EVAL]
                  [--two_gene_cov_seq TWO_GENE_COV_SEQ] [--two_gene_cov_hmm TWO_GENE_COV_HMM] [--single_gene_eval SINGLE_GENE_EVAL] [--single_gene_cov_seq SINGLE_GENE_COV_SEQ]
                  [--single_cov_hmm SINGLE_COV_HMM] [--vf_eval VF_EVAL] [--vf_cov_hmm VF_COV_HMM]
                  input output

positional arguments:
  input                 Input fasta file
  output                Prefix for output directory

optional arguments:
  -h, --help            show this help message and exit
  -t THREADS, --threads THREADS
                        Number of parallel processes [4].
  --prodigal {single,meta}
                        Which mode to run prodigal in [single].
  --aa                  Input is a protein fasta. Has to be in prodigal format.
  --skip_check          Skip check of input.
  --keep_prodigal       Keep prodigal output.
  --log_lvl {DEBUG,INFO,WARNING,ERROR}
                        Logging level [INFO].
  --redo_typing         Redo the typing. Skip prodigal and HMMER and load the hmmer.tab from the output dir.

data arguments:
  --scores SCORES       Path to CasScoring table.
  --hmms HMMS           Path to directory with HMM profiles.

threshold arguments:
  --dist DIST           Max allowed distance between genes in operon [3].
  --overall_eval OVERALL_EVAL
                        Overall E-value threshold [0.001].
  --overall_cov_seq OVERALL_COV_SEQ
                        Overall sequence coverage threshold [0.5].
  --overall_cov_hmm OVERALL_COV_HMM
                        Overall HMM coverage threshold [0.5].
  --two_gene_eval TWO_GENE_EVAL
                        Two-gene operon E-value threshold [1e-05].
  --two_gene_cov_seq TWO_GENE_COV_SEQ
                        Two-gene operon sequence coverage threshold [0.8].
  --two_gene_cov_hmm TWO_GENE_COV_HMM
                        Two-gene operon HMM coverage threshold [0.8].
  --single_gene_eval SINGLE_GENE_EVAL
                        Lonely gene E-value threshold [1e-10].
  --single_gene_cov_seq SINGLE_GENE_COV_SEQ
                        Lonely gene sequence coverage threshold [0.9].
  --single_cov_hmm SINGLE_COV_HMM
                        Lonely gene HMM coverage threshold [0.9].
  --vf_eval VF_EVAL     V-F Cas12 specific E-value threshold [1e-75].
  --vf_cov_hmm VF_COV_HMM
                        V-F Cas12 specific HMM coverage threshold [0.97].
