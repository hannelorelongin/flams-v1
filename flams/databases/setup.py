from pathlib import Path
import subprocess
import os
from flams.databases import cplmv4
from flams.utils import get_data_dir
from dataclasses import dataclass


@dataclass
class ModificationDatabase:
    module: any
    descriptor: str


@dataclass
class ModificationType:
    type: str
    version: float
    dbs: list[ModificationDatabase]


# Here we store a dict of modifications that can be queried for.
# Each modification has a dbs and version attribute.
# Dbs is a list of tuples (module, label) where module is used to get a FASTA of the modifications using label
# TODO Add more modifications from CPLM
MODIFICATIONS = {
    "acetylation": ModificationType(
        "acetylation", 1.0, [ModificationDatabase(cplmv4, "Acetylation")]
    ),
    "lactylation": ModificationType(
        "lactylation", 1.0, [ModificationDatabase(cplmv4, "Lactylation")]
    ),
    "formylation": ModificationType(
        "formulation", 1.0, [ModificationDatabase(cplmv4, "Formylation")]
    ),
    "succinylation": ModificationType(
        "succinylation", 1.0, [ModificationDatabase(cplmv4, "Succinylation")]
    ),
    "hmgylation": ModificationType(
        "hmglyation", 1.0, [ModificationDatabase(cplmv4, "HMGylation")]
    ),
}


# modifications: is a list of modification names (strings) from user input
def update_db_for_modifications(list_of_mods_to_check: list[str]):
    for m in list_of_mods_to_check:
        _generate_blastdb_if_not_up_to_date(MODIFICATIONS[m])


# Check if data dir contains a BLASTDB with name {modification}-{version}
# If not, download FASTAs and generate BLASTDB.
def _generate_blastdb_if_not_up_to_date(modification: ModificationType):
    BLASTDB_PATH = get_blastdb_path_for_modification(
        modification.type, modification.version
    )

    # If an up-to-date BLASTDB for the given modification already exists, do nothing.
    if Path(f"{BLASTDB_PATH}.pdb").exists():
        return

    fasta_location = f"{get_data_dir()}/{modification.type}.fasta"

    if os.path.exists(fasta_location):
        os.remove(fasta_location)

    _get_fasta_from_dbs(modification, fasta_location)

    # Generate local BLASTDB from FASTA in fasta_location
    _generate_blastdb(fasta_location, BLASTDB_PATH)


def _get_fasta_from_dbs(modification: ModificationType, fasta_location):
    for db in modification.dbs:
        db.module.get_fasta(db.descriptor, fasta_location)


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


# Gets the path to where the local BLASTDB for a given modification should be located.
def get_blastdb_path_for_modification(modification: str, version=None):
    # If version was not specified, get the current
    if not version:
        version = MODIFICATIONS[modification].version

    return f"{get_data_dir()}/{modification}-{version}"
