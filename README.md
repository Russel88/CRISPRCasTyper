# CasPredict

Detect CRISPR-Cas genes, group them into operons, and predict their subtype 

## Installation
### Requirements
##### Prodigal, hmmer, and some python3 modules
It is advised to use [miniconda](https://docs.conda.io/en/latest/miniconda.html) or [anaconda](https://www.anaconda.com/) to manage dependencies.
However, if you have the dependencies in your PATH it will still work.
```sh
conda create -n caspredict python=3.8 prodigal hmmer biopython pandas scipy multiprocess
```

grep and sed should also be in your PATH. If you don't already have these, include them in your conda environment:
```sh
conda create -n caspredict python=3.8 prodigal hmmer biopython pandas scipy multiprocess grep sed
```

### CasPredict
##### Clone git repo
```sh
git clone https://github.com/Russel88/CasPredict.git
```

## How to run
##### Activate environment
```sh
conda activate caspredict
```
##### Run with a nucleotide fasta as input
```sh
./CasPredict.py genome.fa my_output
```
##### Use multiple threads
```sh
./CasPredict.py genome.fa my_output -t 20
```
##### Check the different options
```sh
./CasPredict.py -h

usage: CasPredict.py [-h] [-t THREADS] [--prodigal PRODIGAL] [--aa] [--check_input CHECK_INPUT] [--keep_prodigal KEEP_PRODIGAL] [--log_lvl LOG_LVL] [--scores SCORES] [--hmms HMMS] [--dist DIST]
                     [--overall_eval OVERALL_EVAL] [--overall_cov_seq OVERALL_COV_SEQ] [--overall_cov_hmm OVERALL_COV_HMM] [--two_gene_eval TWO_GENE_EVAL] [--two_gene_cov_seq TWO_GENE_COV_SEQ]
                     [--two_gene_cov_hmm TWO_GENE_COV_HMM] [--single_gene_eval SINGLE_GENE_EVAL] [--single_gene_cov_seq SINGLE_GENE_COV_SEQ] [--single_cov_hmm SINGLE_COV_HMM] [--vf_eval VF_EVAL]
                     [--vf_cov_hmm VF_COV_HMM]
                     input output

positional arguments:
  input                 Input fasta file
  output                Prefix for output directory

optional arguments:
  -h, --help            show this help message and exit
  -t THREADS, --threads THREADS
                        Number of parallel processes. Default 4
  --prodigal PRODIGAL   Which mode to run prodigal in. Default single
  --aa                  Input is a protein fasta. Has to be in prodigal format
  --check_input CHECK_INPUT
                        Should the input be checked. Default True
  --keep_prodigal KEEP_PRODIGAL
                        Keep prodigal output. Default False
  --log_lvl LOG_LVL     Logging level. Default 20
  --redo_typing         Redo the typing. Skip prodigal and HMMER and load the hmmer.tab from the output dir

data arguments:
  --scores SCORES       Path to CasScoring table. Default same dir as CasPredict script
  --hmms HMMS           Path to directory with HMM profiles. Default same dir as CasPredict script

threshold arguments:
  --dist DIST           Max allowed distance between genes in operon. Default 3
  --overall_eval OVERALL_EVAL
                        Overall E-value threshold. Defalt 0.001
  --overall_cov_seq OVERALL_COV_SEQ
                        Overall sequence coverage threshold. Default 0.5
  --overall_cov_hmm OVERALL_COV_HMM
                        Overall HMM coverage threshold. Default 0.5
  --two_gene_eval TWO_GENE_EVAL
                        Two-gene operon E-value threshold. Default 1e-5
  --two_gene_cov_seq TWO_GENE_COV_SEQ
                        Two-gene operon sequence coverage threshold. Default 0.8
  --two_gene_cov_hmm TWO_GENE_COV_HMM
                        Two-gene operon HMM coverage threshold. Default 0.8
  --single_gene_eval SINGLE_GENE_EVAL
                        Lonely gene E-value threshold. Default 1e-10
  --single_gene_cov_seq SINGLE_GENE_COV_SEQ
                        Lonely gene sequence coverage threshold. Default 0.9
  --single_cov_hmm SINGLE_COV_HMM
                        Lonely gene HMM coverage threshold. Default 0.9
  --vf_eval VF_EVAL     V-F Cas12 specific E-value threshold. Default 1e-75
  --vf_cov_hmm VF_COV_HMM
                        V-F Cas12 specific HMM coverage threshold. Default 0.97
```
