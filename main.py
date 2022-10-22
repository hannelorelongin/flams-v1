import csv
import argparse
from pathlib import Path

from Bio import SeqIO
from Bio.Blast import NCBIWWW
from flams.display import display_result
from flams.run_blast import run_blast


def read_fasta(file: Path):
    record = SeqIO.read(file, "fasta")
    return record

def main(args):
    protein_seq = read_fasta(Path(args.input[0]))
    
    result = run_blast(protein_seq, *args)
    display_result(result)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Find Lysine Acetylation Modification Sites.')
    parser.add_argument("input", nargs=1, type=str, help="Path to input .fasta file.")
    parser.add_argument("-o", "--output", nargs=1, type=str, default="out.tsv", help="Path to output .tsv file.")
    args = parser.parse_args()

    main(args)