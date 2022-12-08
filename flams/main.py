from pathlib import Path
import sys

from flams.display import display_result
from flams.input import parse_args
from flams.run_blast import run_blast
from flams.databases.setup import update_db_for_modifications


DATA_PATH = Path(__file__).parents[1] / "data"


def main(args, protein_file):
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
