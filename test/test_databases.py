from pathlib import Path
import unittest
from unittest.mock import patch
import responses
import glob
from flams.databases.setup import get_blastdb_name_for_modification, update_db_for_modifications
import flams.databases.cplmv4 as cplmv4
from Bio import SeqIO

TEST_OUTPUT_PATH = Path(__file__).parent / "testfiles"
TEST_OUTPUT_FASTA = TEST_OUTPUT_PATH / "hmgylation.fasta.tmp"

TEST_INPUT_HMGYLATION = TEST_OUTPUT_PATH / "HMGylation.zip"


class DatabaseTestCase(unittest.TestCase):
    def tearDown(self) -> None:
        for file in TEST_OUTPUT_PATH.glob("*.tmp"):
            file.unlink()

        blastdb_name = get_blastdb_name_for_modification(
            "hmgylation-unittest", 0.1
        )
        for file in glob.glob(f"{blastdb_name}.*"):
            Path(file).unlink()

    # This tests get_fasta and _convert_plm_to_fasta of the CPLMv4 module
    # Executed by test_fasta_download_and_blastdb_generation
    def test_download_cplmv4_unzip_and_write(self):
        self.assertFalse(TEST_OUTPUT_FASTA.exists())

        with responses.RequestsMock() as rsps:
            with open(TEST_INPUT_HMGYLATION, "rb") as file:
                file_content = file.read()
                rsps.add(
                    responses.GET,
                    cplmv4.URL.format("HMGylation"),
                    body=file_content,
                    status=200,
                    content_type="application/zip",
                    headers={"content-length": str(len(file_content))},
                )
            cplmv4.get_fasta("HMGylation", TEST_OUTPUT_FASTA)

        self.assertTrue(TEST_OUTPUT_FASTA.exists())

        # Validate file is in FASTA format
        with open(TEST_OUTPUT_FASTA, "r") as handle:
            fasta = SeqIO.parse(handle, "fasta")
            # Assert that the HMGylation DB contains 126 entries.
            self.assertTrue(len(list(fasta)) == 126)

    @patch("flams.databases.setup._generate_blastdb", return_value=None)
    @patch("flams.databases.setup._get_fasta_from_dbs", return_value=None)
    @patch("flams.databases.setup.get_data_dir", return_value=TEST_OUTPUT_PATH)
    def test_dont_download_fasta_if_already_exists(self, get_data_dir, get_fasta, generate_blastdb):
        self.assertTrue((TEST_OUTPUT_PATH / "hmgylation-1.0.fasta").exists())

        update_db_for_modifications(["hmgylation"])

        get_fasta.assert_not_called()


