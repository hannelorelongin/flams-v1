#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: kasgel, hannelorelongin, annkamsk
"""

import argparse
import logging
import sys

from pathlib import Path
from typing import Tuple
from Bio import SeqIO
from databases import setup as db_setup
import requests

""" setup
This script deals with parsing the input and checking the validity of all provided arguments.
"""

DATA_PATH = Path(__file__).parent / "data"

logging.basicConfig(
        level = logging.INFO,
        format = '[%(asctime)s] %(levelname)s: %(message)s',
        datefmt = '%d/%m %H:%M:%S',
        force=True
    )


def parse_args(sys_args) -> Tuple[argparse.Namespace, Path]:
    """
    This function parses all provided arguments.

    Parameters
    ----------
    sys_args:
        Arguments passed to FLAMS

    """
    parser = create_args_parser()
    args = parser.parse_args(sys_args)
    protein_file = validate_input(args, parser)
    return args, protein_file


def create_args_parser():
    """
    This function creates an argument parser.

    """
    parser = argparse.ArgumentParser(
        description="Find Lysine Acylation & other Modification Sites."
    )

    # query proteins
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        "--in",
        dest="input",
        type=Path,
        help="Path to input .fasta file.",
        metavar="inputFilePath"
    )
    group.add_argument(
        "--id",
        type=str,
        help="UniProt ID of input protein.",
        metavar="UniProtID"
    )

    # position
    parser.add_argument(
        "-p",
        "--pos",
        required=True,
        type=int,
        help="Position in input protein that will be searched for conserved modifications.",
        metavar="position"
    )

    parser.add_argument(
        "-m",
        "--modification",
        nargs="+",
        default=["Acylations"],
        help="List of modifications to search for at the given lysine position. Possible values  is one or a combination (seperated by spaces) of: ubiquitination, sumoylation, pupylation, neddylation, acetylation, succinylation, crotonylation, malonylation, 2-hydroxyisobutyrylation, beta-hydroxybutyrylation, butyrylation, propionylation, glutarylation, lactylation, formylation, benzoylation, hmgylation, mgcylation, mgylation, methylation, glycation, hydroxylation, phosphoglycerylation, carboxymethylation, lipoylation, carboxylation, dietylphosphorylation, biotinylation, carboxyethylation. We also provide aggregated combinations: 'All','Ubs','Acylations' and'Others', in analogy to the CPLM database. [default: Acylations]",
        metavar="modification"
    )

    parser.add_argument(
        "--range",
        type=int,
        default=0,
        help="Allowed error range for position. [default: 0]",
        metavar="errorRange"
    )

    # BLAST settings
    parser.add_argument(
        "-t",
        "--num_threads",
        type=int,
        default=1,
        help="Number of threads to run BLAST with. [default: 1]",
        metavar="threadsBLAST"
    )

    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=Path("out.tsv"),
        help="Path to output .tsv file. [default: out.tsv]",
        metavar="outputFilePath",
    )

    return parser


def validate_input(args, parser) -> Path:
    """
    This function checks whether all arguments pass the checks.

    Parameters
    ----------
    args:
        Arguments passed to flams
    parser:
        Argument parser

    """
    check_files_valid(args, parser)

    protein_file = get_protein_file(args, parser)

    check_lysine(protein_file, args.pos, parser)
    check_modifications(args, parser)
    return protein_file


def check_files_valid(args, parser):
    """
    This function checks whether (i) the provided input FASTA file exists and is a valid FASTA file and (ii) if the provided output file is a file, not a directory.

    Parameters
    ----------
    args:
        Arguments passed to flams
    parser:
        Argument parser

    """
    if args.input:
        if not args.input.exists():
            logging.error(f"Input file {args.input} does not exist. Please provide the correct path to the input file. Exiting FLAMS...")
            sys.exit()

        if not is_valid_fasta_file(args.input):
            logging.error(f"Input file {args.input} is not a valid FASTA file. Exiting FLAMS...")
            sys.exit()

    if args.output and args.output.is_dir():
        logging.error(f"Provided output: {args.output} is a directory name, not a file name. Please provide an output filename instead. Exiting FLAMS...")
        sys.exit()


def is_valid_fasta_file(path: Path):
    """
    This function checks whether the provided input FASTA file is a valid FASTA file.

    Parameters
    ----------
    path: Path
        Path to input FASTA file containing info on query protein.

    """
    try:
        SeqIO.read(path, "fasta")
        return True
    except Exception:
        return False


def get_protein_file(args, parser) -> Path:
    """
    This function retrieves the protein input file, by either:
    - returning the path to the user provided input file through option --input directly
    - downloading the FASTA file from UniProt, based on the user provided UniProt ID, then returning the path to the downloaded protein fasta

    Parameters
    ----------
    args:
        Arguments passed to flams
    parser:
        Argument parser

    """
    if args.input:
        return args.input

    try:
        return retrieve_protein_from_uniprot(args.id)
    except requests.HTTPError:
        logging.error("Non-existing UniProt ID. Please provide a valid UniProt ID. Exiting FLAMS...")
        sys.exit()


def retrieve_protein_from_uniprot(uniprot_id) -> Path:
    """
    This function downloads the FASTA file from UniProt, based on the provided UniProt ID, then returns the path to the downloaded protein fasta.

    Parameters
    ----------
    uniprot_id: str
        UniProt ID of query protein as provided by user with --id

    """
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta"
    logging.info(f"Retrieving FASTA file for Uniprot ID {uniprot_id} at {url}")
    r = requests.get(url)

    r.raise_for_status()

    filename = DATA_PATH / f"{uniprot_id}.fasta.tmp"
    with filename.open("w+") as f:
        f.write(r.text)

    logging.info(f"Stored FASTA file for Uniprot ID {uniprot_id} at {filename}")
    return filename


def check_lysine(protein_file, pos, parser):
    """
    This function checks whether the user provided position actually points to a lysine.
    If not, it returns an error.

    Parameters
    ----------
    protein_file: fasta
        FASTA file containing query protein
    pos: int
        User provided position in the query protein
    parser:
        Argument parser

    """
    try:
        if not is_position_lysine(pos, protein_file):
            logging.error(
                f"Position {pos} does not point to lysine: {_get_position_display_str(pos, protein_file)} "
                )
            logging.error("Please provide a position that corresponds to a lysine.")
            sys.exit()
    except IndexError as e:
        logging.error(f"{e}. Please provide a lysine position smaller than the size of your protein.")
        sys.exit()

def is_position_lysine(position: int, input: Path) -> bool:
    """
    This function assess whether the user provided position actually points to a lysine in the query protein.

    Parameters
    ----------
    position: int
        User provided position in the query protein
    input: Path
        Path to FASTA file containing query protein

    """
    # user provides position in 1-based indexing system
    position_idx = position - 1
    input_seq = SeqIO.read(input, "fasta").seq
    return input_seq[position_idx] == "K"


def check_modifications(args, parser):
    """
    This function checks whether the user provided modification is part of the collection of modifications that can be queried.
    If not, it returns an error. It also transforms the aggregate modification options Ubs, Acylations, Others, and All
    to their respective collection of modifications. Finally, it removes any duplicate modification types.

    Parameters
    ----------
    args:
        Arguments passed to flams
    parser:
        Argument parser

    """
    if args.modification:
        if ('Ubs' in args.modification) | ('Acylations' in args.modification) | ('Others' in args.modification) | ('All' in args.modification):
            if 'Ubs' in args.modification:
                args.modification.remove('Ubs')
                args.modification.extend(['ubiquitination','sumoylation','pupylation','neddylation'])
            if 'Acylations' in args.modification:
                args.modification.remove('Acylations')
                args.modification.extend(['lactylation','acetylation','succinylation','crotonylation','malonylation',
                'beta-hydroxybutyrylation','benzoylation','propionylation','2-hydroxyisobutyrylation','formylation',
                'hmgylation','mgcylation','mgylation','glutarylation','butyrylation'])
            if 'Others' in args.modification:
                args.modification.remove('Others')
                args.modification.extend(['methylation','hydroxylation','phosphoglycerylation','biotinylation','lipoylation',
                'dietylphosphorylation','glycation','carboxymethylation','carboxyethylation','carboxylation'])
            if 'All' in args.modification:
                args.modification.remove('All')
                args.modification.extend(['ubiquitination','sumoylation','pupylation','neddylation',
                'lactylation','acetylation','succinylation','crotonylation','malonylation',
                'beta-hydroxybutyrylation','benzoylation','propionylation','2-hydroxyisobutyrylation','formylation',
                'hmgylation','mgcylation','mgylation','glutarylation','butyrylation',
                'methylation','hydroxylation','phosphoglycerylation','biotinylation','lipoylation',
                'dietylphosphorylation','glycation','carboxymethylation','carboxyethylation','carboxylation'])
        args.modification = list(set(args.modification))
        for i in args.modification:
            if i not in db_setup.MODIFICATIONS:
                logging.error(f"Invalid modification type {i}. Please choose a modification from the list specified in the docs. ")
                sys.exit()

def _get_position_display_str(position: int, input: Path) -> str:
    """
    This function returns a fragment of the sequence around a chosen position.

    Parameters
    ----------
    position: int
        User provided position in the query protein
    input: Path
        Path to FASTA file containing query protein

    """
    # user provides position in 1-based indexing system
    pos_idx = position - 1
    seq = SeqIO.read(input, "fasta").seq
    lower = max(0, pos_idx - 3)
    upper = min(len(seq), pos_idx + 3)
    prefix = "..." if lower > 0 else ""
    sufix = "..." if upper < len(seq) - 1 else ""
    pos_idx = len(prefix) + (pos_idx - lower)
    seq_row = f"{prefix}{seq[lower:upper]}{sufix}"
    pointer_row = " " * pos_idx + "^"
    return "".join(["\n", seq_row, "\n", pointer_row])
