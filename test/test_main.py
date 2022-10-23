from pathlib import Path
import unittest
from Bio import SeqIO

from flams.main import retrieve_protein_from_uniprot


DATA_PATH = Path(__file__).parents[1] / "data"


class GetProteinTestCase(unittest.TestCase):
    def tearDown(self) -> None:
        for file in DATA_PATH.glob("*.tmp"):
            file.unlink()

    def test_retrieve_from_uniprot(self):
        # given
        protein_id = "P57703"
        expected_filename = DATA_PATH / f"{protein_id}.fasta.tmp"
        self.assertFalse(expected_filename.exists())

        # when
        filename = retrieve_protein_from_uniprot("P57703")

        # then
        self.assertEqual(filename, expected_filename)
        self.assertTrue(expected_filename.exists())
        with filename.open() as f:
            seq = SeqIO.read(f, "fasta")
        self.assertEqual(seq.id, "sp|P57703|METE_PSEAE")
