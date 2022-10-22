from pathlib import Path
import unittest

from Bio.Blast import NCBIXML
from flams.display import display_result


EXAMPLE_BLAST_OUTPUT = Path(__file__).parent / "testfiles/acetylation.out"


class DisplayTestCase(unittest.TestCase):
    def test_display_simple(self):
        with EXAMPLE_BLAST_OUTPUT.open() as f:
            blast_records = list(NCBIXML.parse(f))
            display_result(blast_records)


if __name__ == "__main__":
    unittest.main()
