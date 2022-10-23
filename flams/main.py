from pathlib import Path
import sys

from Bio import SeqIO
from display import display_result
from input import parse_args
from run_blast import run_blast


def read_fasta(file: Path):
    record = SeqIO.read(file, "fasta")
    return record


def main(args):
    #protein_seq = read_fasta(Path(args.input))

    result = run_blast(
        args.input,
        lysine_pos=args.pos,
        lysine_range=args.range,
        num_threads=args.num_threads,
    )
    display_result(result)


if __name__ == "__main__":
    args = parse_args(sys.argv[1:])
    main(args)
