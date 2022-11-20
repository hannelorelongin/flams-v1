import argparse
from pathlib import Path
from Bio import SeqIO
from flams.databases import setup as db_setup


def parse_args(sys_args):
    parser = create_args_parser()
    args = parser.parse_args(sys_args)
    check_files_valid(args, parser)
    check_lysine(args, parser)
    check_modifications(args, parser)
    # TODO(annkamsk) check range
    return args


def create_args_parser():
    parser = argparse.ArgumentParser(
        description="Find Lysine Acetylation Modification Sites."
    )

    # query proteins
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        "--in",
        dest="input",
        type=Path,
        help="Path to input .fasta file.",
    )
    group.add_argument(
        "--id", type=str, help="Uniprot ID of input protein.", metavar="id"
    )

    # position
    parser.add_argument(
        "pos", type=int, help="Position in input where to search for Acetylation"
    )
    parser.add_argument(
        "--range",
        type=int,
        default=0,
        help="Range in position where to search for Acetylation",
    )

    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=Path("out.tsv"),
        help="Path to output .tsv file.",
        metavar="output",
    )

    # TODO add full list of possible values
    parser.add_argument(
        "-m",
        "--modification",
        nargs="+",
        default=["acetylation"],
        help="List of modifications to search for at the given lysine position. Possible values: acetylation, "
        "lactylation, formylation",
    )

    # BLAST settings
    parser.add_argument(
        "-t",
        "--num_threads",
        type=int,
        default=1,
        help="Number of threads to run BLAST with",
    )
    return parser


def check_files_valid(args, parser):
    if not args.input.exists():
        parser.error(f"Input file {args.input} does not exist.")

    if not is_valid_fasta_file(args.input):
        parser.error(f"Input file {args.input} is not a valid fasta file.")

    if args.output and args.output.is_dir():
        parser.error(f"Provided output: {args.output} is a directory.")


def is_valid_fasta_file(path: Path):
    try:
        SeqIO.read(path, "fasta")
        return True
    except Exception:
        return False


def check_lysine(args, parser):
    if not is_position_lysine(args.pos, args.input):
        parser.error(
            f"Position {args.pos} does not point to Lysine: {_get_position_display_str(args.pos, args.input)}"
        )


def is_position_lysine(position: int, input: Path) -> bool:
    # user provides position in 1-based indexing system
    position_idx = position - 1
    input_seq = SeqIO.read(input, "fasta").seq
    return input_seq[position_idx] == "K"


def check_modifications(args, parser):
    if args.modification:
        for i in args.modification:
            if i not in db_setup.MODIFICATIONS:
                parser.error(f"Invalid modification type {i}")


def _get_position_display_str(position: int, input: Path) -> str:
    """
    Display fragment of sequence around chosen position
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
