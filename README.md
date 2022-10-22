# Find Lysine Acetylation Modification Sites (working name)

Requirements:
* python >= 3.10

## Setup
Download the project:
`git clone git@github.com:annkamsk/flams.git`
`cd flams`

Create virtual environment and activate it:
`python -m venv venv/`
`source venv/bin/activate`

Install dependencies:
`pip install -r requirements.txt`

## Run
`python flams/main.py --in {input fasta file} N`
where N is position of lysine

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


# Dev notes

## Creating local BLAST database
Manual: http://biopython.org/DIST/docs/tutorial/Tutorial.html#sec125

`makeblastdb -in data/acetylation.faa -dbtype prot`

## Parsing PLM database format
Go to `read_plm.py`, change `PLM_DATABASE` to path of database, and `OUTPUT` to output path. Then run:
`python read_plm.py`