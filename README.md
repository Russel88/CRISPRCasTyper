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
./CasPredict.py -i genome.fa -o my_output
```
##### Use multiple threads
```sh
./CasPredict.py -i genome.fa -o my_output -t 20
```
##### Check the different options
```sh
./CasPredict.py -h
```
