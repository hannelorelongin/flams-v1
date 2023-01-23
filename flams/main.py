#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: kasgel, Retro212, annkamsk, hannelorelongin
"""

from pathlib import Path
import sys
import shutil

from flams.display import display_result
from flams.input import parse_args
from flams.run_blast import run_blast
from flams.databases.setup import update_db_for_modifications


DATA_PATH = Path(__file__).parents[1] / "data"

def is_available(program):
    """Verify if third party dependency program is available on the PATH"""

    if shutil.which(program) is not None:
        print(program + " installed: OK")
    else:
        print(program + " installed: Not available on the path.. exiting FLAMS.")
        sys.exit()

def main(args, protein_file):

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


if __name__ == "__main__":
    args, input_file = parse_args(sys.argv[1:])
    main(args, input_file)
