# FLAMS: Find Lysine Acylations & other Modification Sites

A bioinformatics tool to analyze the conservation of lysine modifications, by means of a position-based search against the Compendium of Protein Lysine Modifications (CPLM database) v.4. FLAMS is available as command-line tool and as a [web service](https://www.biw.kuleuven.be/m2s/cmpg/research/CSB/tools/flams/).

# Table of contents

1.  [Introduction](#introduction)
2.  [System requirements](#system-requirements)
    1.  [General dependencies](#general-dependencies)
    2.  [Third-party dependencies](#third-party-dependencies)
3.  [Installation](#installation)
4.  [Usage](#usage)
    1. [Example use case](#example-use-case)
5.  [Output](#output)
6.  [Contact](#contact)
7.  [References](#references)
8.  [License](#license)

## Introduction

FLAMS is a bioinformatics tool to analyze the conservation of lysine modifications, by means of a position-based search against the CPLM database v.4 (Zhang, W. et al. Nucleic Acids Research. 2021, 44 (5): 243–250.). FLAMS can be used (i) to quickly verify whether modifications in a specific protein have been reported before, (ii) to assess whether findings in one species might translate to other species, and (iii) to systematically assess the novelty and conservation of all reported lysine modification sites.

The tool takes as input a protein (identifier or sequence) and the position of a lysine. This repository contains the command-line tool `FLAMS`, which obtains an overview of the previously reported lysine modifications matching your query, by using the following scripts:

* *input.py*: processing the user-provided input
* *cplm4.py* and *setup.py*: downloading and preparing the modification-specific databases
* *run_blast.py*: searching your query against the databases of proteins with lysine modifications
* *display.py*: formatting the list of conserved lysine modifications to a tab delimited output file

FLAMS is also available as a web service at https://www.biw.kuleuven.be/m2s/cmpg/research/CSB/tools/flams/ .

## System requirements

Linux 64-bit, Windows and Mac OS supported.

### General dependencies

* Python3 (>=3.10)

### Third-party dependencies

* [BLAST+ v.13](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)

## Installation

The recommended installation for Mac OS and Linux is through conda:

`conda install -c bioconda flams`

It is also possible to install FLAMS through pip (recommended installation for Windows):

`pip install flams`

Please note that the pip install requires users to have BLAST+ installed locally and available in PATH. For more information on how to install BLAST+ on Windows, click [here](https://www.ncbi.nlm.nih.gov/books/NBK52637/) .

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

### Example use case

We provide two example use cases for FLAMS:

With the following command, you search whether the TatA (UniProt ID: A0A916NWA0) acetylation on K66 in *Dehalococcoide mccartyi* strain CBDB1, as described by [Greiner-Haas (2021)](https://doi.org/10.3390/microorganisms9020365), had been previously detected.

`FLAMS --in A0A916NWA0.fa -p 66 -m acetylation -o tatA.tsv`

With the following command, you search whether the *Mycobabcterium smegmatis*' FadD2 (UniProt ID: A0QQ22) K537 is known to carry any modifications of the 'acylations' category, similar to what was reported by [Xu (2020)](https://doi.org/10.1128/mSystems.00424-19).

`FLAMS --id A0QQ22 -p 537 -m Acylations -o FadD2.tsv`

You can find the example input and output data in the folder `test_data`.

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

If you use FLAMS in your work, please cite us.

In addition, FLAMS relies on third-party software & database:

Zhang, W., Tan, X., Lin, S., Gou, Y., Han, C., Zhang, C., Ning, W., Wang, C. & Xue, Y. (2021) "CPLM 4.0: an updated database with rich annotations for protein lysine modifications." Nucleic Acids Research. 44(5):243–250.

Altschul, S.F., Gish, W., Miller, W., Myers, E.W. & Lipman, D.J. (1990) "Basic local alignment search tool." J. Mol. Biol. 215:403-410.

## License

FLAMS is freely available under an MIT license.

Use of the third-party software, libraries or code referred to in the References section above may be governed by separate terms and conditions or license provisions. Your use of the third-party software, libraries or code is subject to any such terms and you should check that you can comply with any applicable restrictions or terms and conditions before use.
