# FLAMS: Find Lysine Acetylation Modification Sites (working name)

A bioinformatics tool to analyze the conservation of lysine modifications,
by means of a position-based search against the CPLM database v.4.

# Table of contents

1.  [Introduction](#introduction)
2.  [System requirements](#system-requirements)
    1.  [General dependencies](#general-dependencies)
    2.  [Third-party tools](#third-party-tools)
    3.  [Installation](#installation)
    4.  [Command line options]()
    5.  [Usage](#usage)
    6.  [Output]()
3.  [Contact]()
4.  [References]()

## Introduction

TBD

## System requirements

Linux 64-bit (and OS X?) supported.

### General dependencies

* Python3 (>=3.10)
* shutil
* BioPython
* TBD

### Third-party dependencies

* [BLAST+ v.13](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)

### Installation

The recommended installation takes place in a dedicated conda environment,
which can be created as follows:

`conda create -n flamsEnv`

In this environment, the correct Python environment can be set up,
and other dependencies can be installed through pip:

`conda activate flamsEnv`

`conda install -c conda-forge python=3.10`

`pip install -r requirements.txt`


### Usage

Download the project:

`git clone git@github.com:annkamsk/flams.git`

`cd flams`

Run the tool:

`python -m flams.main [-h] (--in INPUT | --id id) [--range RANGE] [-o output] [-m MODIFICATION [MODIFICATION ...]] [-t NUM_THREADS] pos`

Required arguments:
* `MODIFICATION` is one of: acetylation, lactylation, formylation, succinylation, hmgylation
* one of:
  * `INPUT` is a path to fasta file with the queried protein (has to contain only 1 protein)
  * `id` is a uniprot ID of the queried protein
* `pos` is a position of lysine in the queried protein sequence where you look for modifications

Optional arguments:
* `RANGE` (default: 0) is an number of positions before and after `pos` to also search for modifications
* `NUM_THREADS` (default: 1) is a BLAST parameter
* `output` (default: out.tsv) is path to a csv file where the result will be saved

Example:

`python -m flams.main --in P57703.fa --range 5 -o results.tsv 308 -m acetylation succinylation`

`python -m flams.main --id P57703 19 -m lactylation`

----------------------------------------------------------------------------------------------

# Development

## Linters
Before commiting the code run:
`black .`
`flake8 flams`

## Tests
To run all tests:
`python -m unittest discover`

To run a specific module:
`python -m unittest test.test_display`

## Push a commit
First, create a new branch:
`git checkout -b <new-branch>`
<new-branch> should be a short (1-3 words) hyphen-separated name vaguely related to what you've been working on (eg. `input-read` etc). Don't stress too much about it.
You'll be moved automatically to that branch.

`git add .`  
`git commit -m "{Short description of change}"`  
`git push`  

To merge the code from the branch to the main branch, you need to create a pull request (can be done through the web interface).

## Creating local BLAST database
Manual: http://biopython.org/DIST/docs/tutorial/Tutorial.html#sec125

`makeblastdb -in data/acetylation.faa -dbtype prot`

## Parsing PLM database format
Go to `read_plm.py`, change `PLM_DATABASE` to path of database, and `OUTPUT` to output path. Then run:
`python read_plm.py`
