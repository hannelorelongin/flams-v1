import unittest
from flams.run_blast import _parse_fasta_title


class DisplayTestCase(unittest.TestCase):
    def test_parse_title(self):
        title = "PLMD-7244|P25665|304 Acetylation [Escherichia coli (strain K12)]"
        spl = _parse_fasta_title(title)
        self.assertTrue(spl[0] == "PLMD-7244")
        self.assertTrue(spl[1] == "P25665")
        self.assertTrue(spl[2] == "304")
        self.assertTrue(spl[3] == "Acetylation")
        self.assertTrue(spl[4] == "[Escherichia coli (strain K12)]")
