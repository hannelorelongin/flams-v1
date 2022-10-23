from pathlib import Path
import subprocess
import os
from flams.databases import cplmv4
from flams.utils import get_data_dir

# Here we store a dict of modifications that can be queried for.
# Each modification has a dbs and version attribute.
# Dbs is a list of tuples (module, label) where module is used to get a FASTA of the modifications using label
# TODO Add more modifications from CPLM
MODIFICATIONS = {
    "acetylation": {
        "version": 1.0,
        "dbs": [(cplmv4, "Acetylation")]
    },
    "lactylation": {
        "version": 1.0,
        "dbs": [(cplmv4, "Lactylation")]
    },
    "formylation": {
        "version": 1.0,
        "dbs": [(cplmv4, "Formylation")]
    },
    "succinylation": {
        "version": 1.0,
        "dbs": [(cplmv4, "Succinylation")]
    },
    "hmgylation": {
        "version": 1.0,
        "dbs": [(cplmv4, "HMGylation")]
    }
}


def update_db_for_modifications(modifications):
    for m in modifications:
        _check_if_db_up_to_date(m)


# Check if data dir contains a BLASTDB with name {modification}-{version}
# If not, download FASTAs and generate BLASTDB.
def _check_if_db_up_to_date(modification):
    BLASTDB_PATH = get_blastdb_for_modification(modification)
    if not Path(f'{BLASTDB_PATH}.pdb').exists():
        fasta_location = f'{get_data_dir()}/{modification}.fasta'

        if os.path.exists(fasta_location): os.remove(fasta_location)

        _get_fasta_from_dbs(modification, fasta_location)

        # Generate local BLASTDB from FASTA in fasta_location
        _generate_blastdb(fasta_location, BLASTDB_PATH)


def _get_fasta_from_dbs(modification, fasta_location):
    for db, label in MODIFICATIONS[modification]["dbs"]:
        db.get_fasta(label, fasta_location)


def _generate_blastdb(fasta_in, blastdb_out):
    # Generate local BLASTDB
    try:
        subprocess.call(
            f"makeblastdb -in {fasta_in} -dbtype prot -input_type fasta -parse_seqids -out {blastdb_out}",
            shell=True,
        )
    except FileNotFoundError as e:
        raise FileNotFoundError(
            "You need to install BLAST and include it in system PATH."
        ) from e


def get_blastdb_for_modification(modification):
    version = MODIFICATIONS[modification]["version"]
    return f'{get_data_dir()}/{modification}-{version}'
