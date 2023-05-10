#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: kasgel, Retro212, annkamsk, hannelorelongin
"""

from pathlib import Path
import sys
import shutil
import logging

from .display import display_result
from .input import parse_args
from .run_blast import run_blast
from .databases.setup import update_db_for_modifications

""" FLAMS
FLAMS is a tool that looks for conservation of modifications on lysine residues
by first looking for similar proteins in the Compendium of Protein Lysine Modification Sites (CPLM v.4, Zhang, W. et al. Nucleic Acids Research. 2021, 44 (5): 243–250.),
and then extracting those proteins that contain a modified lysine at the queried position.
The aim of this analysis is to easily identify conserved lysine modifications,
to aid in identifying functional lysine modification sites and the comparison of the results of PTM identification studies across species.

The tool takes as input a protein sequence and the position of a lysine.
"""

def is_available(program):
    """
    This function verifies the installation of third-party dependencies and prints out the result to users.

    Parameters
    ----------
    program: program
        Third-party dependency program.

    """

    if shutil.which(program) is not None:
        logging.info("Checking third-party depencies. Installation of " + program + " : OK.")
    else:
        logging.error("Checking third-party depencies. Installation of " + program + " failed verification: it is not available on the path.. exiting FLAMS.")
        sys.exit()

def main():
    """ Main function of FLAMS
    """


    args, protein_file = parse_args(sys.argv[1:])

    logging.basicConfig(
        level = logging.INFO,
        format = '[%(asctime)s] %(levelname)s: %(message)s',
        datefmt = '%d/%m %H:%M:%S',
        force=True
    )

    is_available('blastp')

    update_db_for_modifications(args.modification)

    # Save absolute path to output file, because run_blast will change working directory.
    output_file = args.output.absolute()

    result = run_blast(
        input=protein_file,
        modifications=args.modification,
        lysine_pos=args.pos,
        lysine_range=args.range,
        num_threads=args.num_threads,
    )

    display_result(output_filename=output_file, blast_records=result)

    logging.info("Succesfully ran FLAMS! You can find your results at: " + str(output_file))


if __name__ == "__main__":
    main()
