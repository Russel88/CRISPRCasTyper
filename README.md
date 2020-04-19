# CRISPRCasTyper

Detect CRISPR-Cas genes and arrays, and predict the subtype based on both Cas genes and CRISPR repeat sequence.

[CRISPRCasTyper and RepeatType are also available through a webserver](http://crisprcastyper.crispr.dk)

This software finds Cas genes with a large suite of HMMs, then groups these HMMs into operons, and predicts the subtype of the operons based on a scoring scheme.
Furthermore, it finds CRISPR arrays with [minced](https://github.com/ctSkennerton/minced), and using a kmer-based machine learning approach (extreme gradient boosting trees) it predicts the subtype of the CRISPR arrays based on the consensus repeat. 
It then connects the Cas operons and CRISPR arrays, producing as output:
* CRISPR-Cas loci, with consensus subtype prediction based on both Cas genes (mostly) and CRISPR consensus repeats
* Orphan Cas operons, and their predicted subtype
* Orphan CRISPR arrays, and their predicted associated subtype

#### It includes the following subtypes:
* All the ones in the most recent Nature Reviews Microbiology (Makarova et al. 2020): [Evolutionary classification of CRISPR–Cas systems: a burst of class 2 and derived variants](https://doi.org/10.1038/s41579-019-0299-x)
* Updated type IV subtypes and variants based on: [Type IV CRISPR–Cas systems are highly diverse and involved in competition between plasmids](https://doi.org/10.1093/nar/gkz1197)
* Type V-K: [RNA-guided DNA insertion with CRISPR-associated transposases](https://doi.org/10.1126/science.aax9181)
* Transposon associated type I-F: [Transposon-encoded CRISPR–Cas systems direct RNA-guided DNA integration](https://doi.org/10.1038/s41586-019-1323-z)

#### It can automatically draw gene maps of CRISPR-Cas systems and orphan Cas operons and CRISPR arrays
<img src='img/plot.png' align="left" height="200" />

#### Citation
Coming soon...

# Table of contents
1. [Quick start](#quick)
2. [Installation](#install)
3. [CRISPRCasTyper - How to](#cctyperhow)
    * [Plotting](#plot)
4. [RepeatType - How to](#repeattype)
4. [RepeatType - Train](#repeattrain)

## Quick start <a name="quick"></a>

```sh
conda create -n cctyper -c conda-forge -c bioconda -c russel88 cctyper
conda activate cctyper
cctyper my.fasta my_output
```

## Installation <a name="install"></a>
CRISPRCasTyper can be installed either through conda or pip.

It is advised to use conda, since this installs CRISPRCasTyper and all dependencies, and downloads with database in one go.

### Conda
Use [miniconda](https://docs.conda.io/en/latest/miniconda.html) or [anaconda](https://www.anaconda.com/) to install.

Create the environment with CRISPRCasTyper and all dependencies and database
```sh
conda create -n cctyper -c conda-forge -c bioconda -c russel88 cctyper
```

### pip
If you have the dependencies (Python >= 3.8, HMMER >= 3.2, Prodigal >= 2.6, grep, sed) in your PATH you can install with pip

```sh
python -m pip install cctyper
```

#### When installing with pip, you need to download the database manually: 
```sh
# Download and unpack
svn checkout https://github.com/Russel88/CRISPRCasTyper/trunk/data
tar -xvzf data/Profiles.tar.gz
mv Profiles/ data/
rm data/Profiles.tar.gz

# Tell CRISPRCasTyper where the data is:
# either by setting an environment variable (has to done for each terminal session, or added to .bashrc):
export CCTYPER_DB="/path/to/data/"
# or by using the --db argument each time you run CRISPRCasTyper:
cctyper input.fa output --db /path/to/data/
```

## CRISPRCasTyper - How to <a name="cctyperhow"></a>
CRISPRCasTyper takes as input a nucleotide fasta, and produces outputs with CRISPR-Cas predictions

#### Activate environment
```sh
conda activate cctyper
```

#### Run with a nucleotide fasta as input
```sh
cctyper genome.fa my_output
```

#### Use multiple threads
```sh
cctyper genome.fa my_output -t 20
```

#### Check the different options
```sh
cctyper -h
```

#### Output <a name="cctyperout"></a>
* **CRISPR_Cas.tab:**           CRISPR_Cas loci, with consensus subtype prediction
    * Contig: Sequence accession
    * Operon: Operon ID (Sequence accession @ NUMBER)
    * Operon_Pos: [Start, End] of operon
    * Prediction: Consenus prediction based on both Cas operon and CRISPR arrays
    * CRISPRs: CRISPRs adjacent to Cas operon
    * Distances: Distances to CRISPRs from Cas operon
    * Prediction_Cas: Subtype prediction based on Cas operon
    * Prediction_CRISPRs: Subtype prediction of CRISPRs based on CRISPR repeat sequences
* **cas_operons.tab:**          All certain Cas operons
    * Contig: Sequence accession
    * Operon: Operon ID (Sequence accession @ NUMBER)
    * Start: Start of operon
    * End: End of operon
    * Prediction: Subtype prediction
    * Best_type: Subtype with the highest score. If the score is high then Prediction = Best_type
    * Best_score: Score of the highest scoring subtype
    * Genes: List of Cas genes
    * Positions: List of Gene IDs for the genes
    * E-values: List of E-values for the genes
    * CoverageSeq: List of sequence coverages for the genes
    * CoverageHMM: List of HMM coverages for the genes
* **crisprs_all.tab:**          All CRISPR arrays
    * Contig: Sequence accession
    * CRISPR: CRISPR ID (Sequence accession _ NUMBER)
    * Start: Start of CRISPR
    * End: End of CRISPR
    * Consensus_repeat: Consensus repeat sequence
    * N_repeats: Number of repeats
    * Prediction: Prediction of the associated subtype based on the repeat sequence
    * Subtype: Subtype with highest prediction probability. Prediction = Subtype if Subtype_probability is high
    * Subtype_probability: Probability of subtype prediction
* **crisprs_orphan.tab:**       Orphan CRISPRs (those not in CRISPR_Cas.tab)
    * Same columns as crisprs_all.tab
* **cas_operons_orphan.tab:**   Orphan Cas operons (those not in CRISPR_Cas.tab)
    * Same columns as cas_operons.tab
* **CRISPR_Cas_putative.tab:**  Putative CRISPR_Cas loci, often lonely Cas genes next to a CRISPR array
    * Same columns as CRISPR_Cas.tab
* **cas_operons_putative.tab:** Putative Cas operons, mostly false positives, but also some ambiguous and partial systems
    * Same columns as cas_operons.tab
* **spacers/*.fa:**             Fasta files with all spacer sequences
* **hmmer.tab:**                All HMM vs. ORF matches, unfiltered results
    * Hmm: HMM name
    * ORF: ORF name (Sequence accession _ Gene ID)
    * tlen: ORF length
    * qlen: HMM length
    * Eval: E-value of alignment
    * score: Alignment score
    * start: ORF start
    * end: ORF end
    * Acc: Sequence accession
    * Pos: Gene ID
    * Cov_seq: Sequence coverage
    * Cov_hmm: HMM coverage
    * strand: Leading (1) or lagging (-1) strand
* **genes.tab**                 All genes and their positions
    * Contig: Sequence accession
    * Start: Start of ORF
    * End: End of ORF
    * Strand: Leading (1) or lagging (-1) strand
    * Pos: Gene ID
* **arguments.tab:**            File with arguments given to CRISPRCasTyper
* **hmmer.log**                 Error messages from HMMER (only produced if any errors were encountered)

##### If run with `--keep_tmp` the following is also produced
* **prodigal.log**              Log from prodigal
* **proteins.faa**              Protein sequences
* **hmmer/*.tab**               Alignment output from HMMER for each Cas HMM
* **minced.out:**               CRISPR array output from minced

#### Notes on output
Files are only created if there is any data. For example, the CRISPR_Cas.tab file is only created if there are any CRISPR-Cas loci. 

### Plotting <a name="plot"></a>
CRISPRCasTyper will automatically plot a map of the CRISPR-Cas loci, orphan Cas operons, and orphan CRISPR arrays.

These maps can be expanded (`--expand N`) by adding unknown genes and genes with alignment scores below the thresholds. This can help in identify potentially un-annotated genes in operons. You can generate new plots without having to re-run the entire pipeline by adding `--redo_typing` to the command. This will re-use the mappings and re-type the operons and re-make the plot, based on new thresholds and plot parameters.

The plot below is run with `--expand 5`

* Cas genes are in red.
* Cas genes, with alignment scores below the thresholds, are in dark green
* Unknown genes are in gray (the number matches the genes.tab file)
* Arrays are in blue, with their predicted subtype association based on the consensus repeat sequence.

<img src='img/plot2.png' align="left" height="350" />

## RepeatType - How to <a name="repeattype"></a>
With an input of CRISPR repeats (one per line, in a simple textfile) RepeatType will predict the subtype, based on the kmer composition of the repeat

#### Activate environment
```sh
conda activate cctyper
```

#### Run with a simple textfile, containing only CRISPR repeats (in capital letters), one repeat per line.
```sh
repeatType repeats.txt
```

#### Output <a name="repeattypeout"></a>
The script prints:
* Repeat sequence
* Predicted subtype
* Probability of prediction

#### Notes on output
* Predictions with probabilities below 0.75 are uncertain, and should be taken with a grain of salt.
* The classifier was only trained on the subtypes for which there were enough (>20) repeats. It can therefore only predict subtypes of repeats associated with the following subtypes:
    * I-A, I-B, I-C, I-D, I-E, I-F, I-G
    * II-A, II-B, II-C
    * III-A, III-B, III-C, III-D
    * IV-A1, IV-A2, IV-A3
    * V-A
    * VI-B
* This is the accuracy per subtype (on an unseen test dataset):
    * **I-A**      0.60
    * **I-B**      0.90
    * **I-C**      0.98
    * **I-D**      0.47
    * **I-E**      1.00
    * **I-F**      0.99
    * **I-G**      0.83
    * **II-A**     0.94
    * **II-B**     1.00
    * **II-C**     0.89
    * **III-A**    0.89
    * **III-B**    0.49
    * **III-C**    0.60
    * **III-D**    0.28
    * **IV-A1**    0.79
    * **IV-A2**    0.78
    * **IV-A3**    0.98
    * **V-A**      0.77
    * **VI-B**     1.00

## RepeatType - Train <a name="repeattrain"></a>
You can train the repeat classifier with your own set of subtyped repeats. With a tab-delimeted input where 1. column contains the subtypes and 2. column contains the CRISPR repeat sequences, RepeatTrain will train a CRISPR repeat classifier that is directly usable for both RepeatType and CRISPRCasTyper.

#### Train
```sh
repeatTrain typed_repeats.tab my_classifier
```

#### Use new model in RepeatType
```sh
repeatType repeats.txt --db my_classifier
```

#### Use new model in CRISPRCasTyper
Save the original database files:
```sh
mv ${CCTYPER_DB}/type_dict.tab ${CCTYPER_DB}/type_dict_orig.tab
mv ${CCTYPER_DB}/xgb_repeats.model ${CCTYPER_DB}/xgb_repeats_orig.model
```

Move the new model into the database folder
```sh
mv my_classifier/* ${CCTYPER_DB}/
```

##### CRISPRCasTyper and RepeatType will now use the new model for repeat prediction!


