import argparse
from pathlib import Path

from Bio import SeqIO
from display import display_result
from run_blast import run_blast


def read_fasta(file: Path):
    record = SeqIO.read(file, "fasta")
    return record


def main(args):
    protein_seq = read_fasta(Path(args.input))

    result = run_blast(
        protein_seq,
        input=args.input,
        lysine_pos=args.lysine_pos,
        lysine_range=args.lysine_range,
        num_threads=args.num_threads,
    )
    display_result(result)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Find Lysine Acetylation Modification Sites."
    )
    parser.add_argument("input", type=str, help="Path to input .fasta file.")
    parser.add_argument(
        "lysine_pos", type=str, help="Position in input where to search for Acetylation"
    )
    parser.add_argument(
        "lysine_range",
        type=str,
        help="Range in position where to search for Acetylation",
    )

    parser.add_argument(
        "-o", "--output", type=str, default="out.tsv", help="Path to output .tsv file."
    )

    # BLAST settings
    parser.add_argument(
        "-t",
        "--num_threads",
        type=int,
        default=1,
        help="Number of threads to run BLAST with",
    )
    args = parser.parse_args()
    main(args)
