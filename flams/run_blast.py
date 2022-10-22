from pathlib import Path
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Blast import NCBIXML
import subprocess

FASTADB_LOCATION = "data/acetylation.faa"
BLASTDB_LOCATION = "data/blast.db"
BLAST_OUT = "data/temp.xml"


def run_blast(
    seq, input, lysine_pos, lysine_range, evalue=0.01, num_threads=1, **kwargs
):
    # First, we need to create a local BLAST DB if it does not exist.
    # we should move this to some user cache directory using e.g. https://pypi.org/project/appdirs/ later

    if not Path(FASTADB_LOCATION).is_file():
        raise FileNotFoundError("Acetylations FASTA file not found.")

    if not Path(f"{BLASTDB_LOCATION}.pdb").is_file():
        try:
            subprocess.call(
                f"makeblastdb -in data/acetylation.faa -dbtype prot -input_type fasta -parse_seqids -out {BLASTDB_LOCATION}",
                shell=True,
            )
        except FileNotFoundError as e:
            raise FileNotFoundError(
                "You need to install BLAST and include it in system PATH."
            ) from e

    # Run BLAST
    blast_exec = NcbiblastpCommandline(
        query=input,
        db=BLASTDB_LOCATION,
        evalue=evalue,
        outfmt=5,
        out=BLAST_OUT,
        num_threads=num_threads,
    )
    blast_exec()

    with open(BLAST_OUT) as handle:
        blast_records = list(NCBIXML.parse(handle))

    return blast_records
