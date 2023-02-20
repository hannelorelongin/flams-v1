# FLAMS: Find Lysine Acylation & other Modification Sites

A command-line tool to analyze the conservation of lysine modifications, by means of a position-based search against the CPLM database v.4.

# Table of contents

1.  [Introduction](#introduction)
2.  [System requirements](#system-requirements)
    1.  [General dependencies](#general-dependencies)
    2.  [Third-party tools](#third-party-tools)
3.  [Installation](#installation)
4.  [Usage](#usage)
5.  [Output](#output)
6.  [Contact](#contact)
7.  [References](#references)

## Introduction

FLAMS is a tool that looks for conservation of modifications on lysine residues by first looking for similar proteins in the Compendium of Protein Lysine Modification Sites (CPLM v.4, Zhang, W. *et al.* Nucleic Acids Research. 2021, 44 (5): 243–250.), and then extracting those proteins that contain a modified lysine at the queried position. The aim of this analysis is to easily identify conserved lysine modifications, to aid in identifying functional lysine modification sites and the comparison of the results of PTM identification studies across species.

The tool takes as input a protein sequence and the position of a lysine. This repository contains a command-line tool `FLAMS` to obtain an overview of the conserved lysine modifications matching your query, by using the following scripts:

* *input.py*: processing the user-provided input
* *cplm4.py* and *setup.py*: downloading and preparing the modification-specific databases
* *run_blast.py*: searching your query against the databases of proteins with lysine modifications
* *display.py*: formatting the list of conserved lysine modifications to a tab delimeted output file

## System requirements

Linux 64-bit supported. Windows users are advised to run the tool through Anaconda.

### General dependencies

* Python3 (>=3.10)
* shutil
* BioPython

### Third-party dependencies

* [BLAST+ v.13](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)

## Installation

The recommended installation takes place in a dedicated conda environment,
which can be created as follows:

`conda create -n flamsEnv`

In this environment, the correct Python environment can be set up,
and other dependencies can be installed through pip:

`conda activate flamsEnv`

`conda install -c conda-forge python=3.10`

`pip install -r requirements.txt`

Please make sure that BLAST is installed locally, and available in the PATH.

## Usage

Download the project:

`git clone git@github.com:annkamsk/flams.git`

`cd flams`

Run the tool:

`python -m flams.main [-h | --help] (--in INPUT | --id id) [--range RANGE] [-o output] [-m MODIFICATION [MODIFICATION ...]] [-t NUM_THREADS] pos`

Required arguments:
* `MODIFICATION` is one or a combination (seperated by spaces) of: ubiquitination, sumoylation, pupylation, neddylation, acetylation, succinylation, crotonylation, malonylation, 2-hydroxyisobutyrylation, beta-hydroxybutyrylation, butyrylation, propionylation, glutarylation, lactylation,  formylation, benzoylation, hmgylation, mgcylation, mgylation, methylation, glycation, hydroxylation, phosphoglycerylation, carboxymethylation, lipoylation, carboxylation, dietylphosphorylation, biotinylation, carboxyethylation
* one of:
  * `INPUT` is the path to a fasta file with the protein you wish to query against (has to contain only 1 protein)
  * `id` is the UniProt ID of the protein you wish to query against
* `pos` is the position of a lysine in the protein, which you want to query against

Optional arguments:
* `RANGE` (default: 0) is an number of positions before and after `pos` to also search for modifications
* `NUM_THREADS` (default: 1) is a BLAST parameter, allows you to speed up the search by multithreading
* `output` (default: out.tsv) is the path to where the result will be saved (tsv file format)

Example:

`python -m flams.main --in P57703.fa --range 5 -o results.tsv 306 -m acetylation succinylation`

`python -m flams.main --id P57703 738 -m acetylation propionylation`

## Output

The output file is a tsv containing one row per modification that matched the query, i.e., aligning (within the user-specified range) to the query lysine, in a protein similar to the query protein. TSV output file contains five columns:
* UniProt ID: UniProt identifier of matched protein
* Species: the textual description of the species of the matched protein
* Modification: the type of modification found in the matched protein
* Lysine location: the location of this matched modification in the matched protein
* Lysine window: the local sequence containing the conserved lysine modification (window of five amino acids before and after)

## Contact

Laboratory of Computational Systems Biology, KU Leuven.

## References

Zhang, W., Tan, X., Lin, S., Gou, Y., Han, C., Zhang, C., Ning, W., Wang, C. & Xue, Y. (2021) "CPLM 4.0: an updated database with rich annotations for protein lysine modifications." Nucleic Acids Research. 44(5):243–250.

Altschul, S.F., Gish, W., Miller, W., Myers, E.W. & Lipman, D.J. (1990) "Basic local alignment search tool." J. Mol. Biol. 215:403-410.
