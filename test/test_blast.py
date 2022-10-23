from pathlib import Path
import unittest
from flams.run_blast import run_blast

TEST_FASTA = Path(__file__).parent / "testfiles/P57703.fasta"

class DisplayTestCase(unittest.TestCase):
    def test_blast_P57703(self):
        out = run_blast(TEST_FASTA, 306, 2)
        self.assertTrue(len(out[0].alignments) == 2)
        # Collect uniprot IDs of matches
        found_ids = set()
        for a in out[0].alignments:
            title = a.title.split("|")
            found_ids.add(title[1])

        self.assertTrue("P25665" in list(found_ids))
        self.assertTrue("P05694" in list(found_ids))
        self.assertTrue(len(found_ids) == 2)

    def test_blast_none(self):
        out = run_blast(TEST_FASTA, 300, 0)
        self.assertTrue(len(out[0].alignments) == 0)
        # Collect uniprot IDs of matches
        found_ids = set()
        for a in out[0].alignments:
            title = a.title.split("|")
            found_ids.add(title[1])
        print(found_ids)
        self.assertTrue(len(found_ids) == 0)