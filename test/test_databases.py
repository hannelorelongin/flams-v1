from pathlib import Path
import unittest
from flams.databases.cplmv4 import get_fasta
from Bio import SeqIO

TEST_OUTPUT_PATH = Path(__file__).parent / "testfiles"
TEST_OUTPUT_FASTA = TEST_OUTPUT_PATH / "hmgylation.fasta.tmp"


class DisplayTestCase(unittest.TestCase):
    def tearDown(self) -> None:
        for file in TEST_OUTPUT_PATH.glob("*.tmp"):
            file.unlink()

    def test_download_hmgylation(self):
        get_fasta("HMGylation", TEST_OUTPUT_FASTA)
        self.assertTrue(TEST_OUTPUT_FASTA.exists())

        # Validate file is in FASTA format
        with open(TEST_OUTPUT_FASTA, "r") as handle:
            fasta = SeqIO.parse(handle, "fasta")
            self.assertTrue(
                any(fasta)
            )  # any(fasta) is false if it was not a fasta file
