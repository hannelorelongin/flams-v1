# FLAMS: Find Lysine Acylations & other Modification Sites

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
* *display.py*: formatting the list of conserved lysine modifications to a tab delimited output file

## System requirements

Linux 64-bit and Mac OS supported. Windows users are advised to run the tool through [Anaconda](https://www.anaconda.com/products/distribution).

### General dependencies

* Python3 (>=3.10)

### Third-party dependencies

* [BLAST+ v.13](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)

## Installation

The recommended installation is through pip:

`pip install flams`

Please make sure that BLAST is installed locally, and available in the PATH.

## Usage

Run the tool:

`FLAMS [-h] (--in inputFilePath | --id UniProtID) -p position [--range errorRange] [-t threadsBLAST] [-o outputFilePath] [-d dataDir] [-m modification [modification ...]]  `

Required arguments:
* one of:
  * `inputFilePath` is the path to a .fasta file with the protein you wish to query against (has to contain only 1 protein)
  * `UniProtID` is the UniProt ID of the protein you wish to query against
* `position` is the position of a lysine in the protein, which you want to query against

Optional arguments:
* `errorRange` is an number of positions before and after `pos` to also search for modifications. [default: 0]
* `threadsBLAST` is a BLAST parameter, allows you to speed up the search by multithreading. [default: 1]
* `outputFilePath` is the path to where the result will be saved (in a .tsv file format). [default: out.tsv]
* `dataDir` is the path to directory where intermediate files (the UniProt sequence files) are stored. [default: $PWD/data]"
* `modification` is one or a combination (seperated by spaces) of: ubiquitination, sumoylation, pupylation, neddylation, acetylation, succinylation, crotonylation, malonylation, 2-hydroxyisobutyrylation, beta-hydroxybutyrylation, butyrylation, propionylation, glutarylation, lactylation,  formylation, benzoylation, hmgylation, mgcylation, mgylation, methylation, glycation, hydroxylation, phosphoglycerylation, carboxymethylation, lipoylation, carboxylation, dietylphosphorylation, biotinylation, carboxyethylation. We also provide aggregated combinations: 'All','Ubs','Acylations' and'Others', in analogy to the CPLM database. [default: Acylations]"

Example:

`FLAMS --in P57703.fa -p 306 --range 5 -o results.tsv -m acetylation succinylation`

`FLAMS --id P57703 -p 738`

For more example use cases, see the Supplementary information of the paper.

## Output

The output file is a .tsv containing one row per modification that matched the query, i.e., aligning (within the user-specified range) to the query lysine, in a protein similar to the query protein. The output file contains five columns:
* UniProt ID: UniProt identifier of matched protein
* Modification: the type of modification found in the matched protein
* Lysine location: the location of this matched modification in the matched protein
* Lysine window: the local sequence containing the conserved lysine modification (window of five amino acids before and after°)
* Species: the textual description of the species of the matched protein

°: window can be smaller than the [-5;+5] window if the sequence alignment ends sooner, which can happen for modified lysines near the start/end of the protein

## Contact

Laboratory of Computational Systems Biology, KU Leuven.

## References

Zhang, W., Tan, X., Lin, S., Gou, Y., Han, C., Zhang, C., Ning, W., Wang, C. & Xue, Y. (2021) "CPLM 4.0: an updated database with rich annotations for protein lysine modifications." Nucleic Acids Research. 44(5):243–250.

Altschul, S.F., Gish, W., Miller, W., Myers, E.W. & Lipman, D.J. (1990) "Basic local alignment search tool." J. Mol. Biol. 215:403-410.
