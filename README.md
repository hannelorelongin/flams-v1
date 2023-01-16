# Find Lysine Acetylation Modification Sites (working name)

Requirements:
	## tested v.3.10.8 - should test all from 3.10-3.11
* python >= 3.10
	## BLAST needs to be installed as executable from site, not through conda, tested working with BLAST 2.13.0
	## Make sure BLAST is installed and added to path for users --> add a check dependency script in tool
* BLAST: https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/

# Usage

## Setup
Download the project:

`git clone git@github.com:annkamsk/flams.git`

`cd flams`

Create virtual environment and activate it:
	##Also works when you create a conda environment for it
`python -m venv venv/`

`source venv/bin/activate`

Install dependencies:

`pip install -r requirements.txt`

## Run
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
