import argparse
from pathlib import Path
from Bio import SeqIO


def parse_args(sys_args):
    parser = create_args_parser()
    args = parser.parse_args(sys_args)
    check_files_valid(args, parser)
    check_lysine(args, parser)
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

    try:
        SeqIO.read(args.input, "fasta")
    except Exception:
        parser.error(f"Input file {args.input} is not a valid fasta file.")

    if args.output and args.output.is_dir():
        parser.error(f"Provided output: {args.output} is a directory.")


def check_lysine(args, parser):
    # user provides position in 1-based indexing system
    position_idx = args.pos - 1
    input_seq = SeqIO.read(args.input, "fasta").seq

    if input_seq[position_idx] != "K":
        parser.error(
            f"Position {args.pos} does not point to Lysine: {_get_position_display_str(input_seq, position_idx)}"
        )


def _get_position_display_str(seq, pos):
    # Display fragment of sequence around chosen position
    lower = max(0, pos - 3)
    upper = min(len(seq), pos + 3)
    prefix = "..." if lower > 0 else ""
    sufix = "..." if upper < len(seq) - 1 else ""
    pos_idx = len(prefix) + (pos - lower)
    seq_row = f"{prefix}{seq[lower:upper]}{sufix}"
    pointer_row = " " * pos_idx + "^"
    return "".join(["\n", seq_row, "\n", pointer_row])
