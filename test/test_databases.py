from pathlib import Path
import unittest
import glob
import flams.databases.dummy
import flams.databases.setup as db_setup
from flams.databases.cplmv4 import get_fasta
from Bio import SeqIO

TEST_OUTPUT_PATH = Path(__file__).parent / "testfiles"
TEST_OUTPUT_FASTA = TEST_OUTPUT_PATH / "hmgylation.fasta.tmp"

TEST_INPUT_HMGYLATION = TEST_OUTPUT_PATH / "HMGylation.zip"


class DatabaseTestCase(unittest.TestCase):
    def tearDown(self) -> None:
        for file in TEST_OUTPUT_PATH.glob("*.tmp"):
            file.unlink()

        blastdb_path = db_setup.get_blastdb_path_for_modification(
            "hmgylation-unittest", 0.1
        )
        for file in glob.glob(f"{blastdb_path}.*"):
            Path(file).unlink()

    def test_download_cplmv4(self):
        self.assertFalse(TEST_OUTPUT_FASTA.exists())
        get_fasta("HMGylation", TEST_OUTPUT_FASTA)
        self.assertTrue(TEST_OUTPUT_FASTA.exists())

        # Validate file is in FASTA format
        with open(TEST_OUTPUT_FASTA, "r") as handle:
            fasta = SeqIO.parse(handle, "fasta")
            # Assert that the HMGylation DB contains 126 entries.
            self.assertTrue(len(list(fasta)) == 126)

    def test_blastdb_generation(self):
        # Assert blastdb for hmgylation-unittest does not exist.
        blastdb_path = db_setup.get_blastdb_path_for_modification(
            "hmgylation-unittest", 0.1
        )
        self.assertFalse(Path(f"{blastdb_path}.pdb").exists())

        # Generate blastdb for hmgylation-unittest using test zip file.
        db = db_setup.ModificationDatabase(flams.databases.dummy, TEST_INPUT_HMGYLATION)
        mod = db_setup.ModificationType("hmgylation-unittest", 0.1, [db])
        db_setup._generate_blastdb_if_not_up_to_date(mod)

        # Assert blastdb for hmgylation-unittest exists.
        self.assertTrue(Path(f"{blastdb_path}.pdb").exists())
