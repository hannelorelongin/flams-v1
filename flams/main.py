from pathlib import Path
import sys
import requests

from flams.display import display_result
from flams.input import parse_args
from flams.run_blast import run_blast
from flams.databases.setup import update_db_for_modifications


DATA_PATH = Path(__file__).parents[1] / "data"


def retrieve_protein_from_uniprot(uniprot_id) -> Path:
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta"
    r = requests.get(url)

    filename = DATA_PATH / f"{uniprot_id}.fasta.tmp"
    with filename.open("w+") as f:
        f.write(r.text)

    return filename


def get_protein(args) -> Path:
    if args.input:
        return args.input
    return retrieve_protein_from_uniprot(args.id)


def main(args):
    update_db_for_modifications(args.modification)
    protein_file = get_protein(args)

    result = run_blast(
        input=protein_file,
        modifications=args.modification,
        lysine_pos=args.pos,
        lysine_range=args.range,
        num_threads=args.num_threads,
    )

    display_result(result)


if __name__ == "__main__":
    args = parse_args(sys.argv[1:])
    main(args)
