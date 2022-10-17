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
`python main.py {input fasta file} -o {output path}`

## Push a commit
First, create a new branch:
`git checkout -b <new-branch>`
<new-branch> should be a short (1-3 words) hyphen-separated name vaguely related to what you've been working on (eg. `input-read` etc). Don't stress too much about it. 
You'll be moved automatically to that branch. 

`git add .`
`git commit -m "{Short description of change}"`
`git push`
