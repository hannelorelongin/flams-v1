from pathlib import Path
import unittest

from Bio.Blast import NCBIXML
from flams.display import display_result


EXAMPLE_BLAST_OUTPUT = Path("test/testfiles/acetylation.out")


def _get_example_blast_output():
    with EXAMPLE_BLAST_OUTPUT.open() as f:
        blast_records = list(NCBIXML.parse(f))
    return blast_records


class DisplayTestCase(unittest.TestCase):
    def test_display_simple(self):
        blast_records = _get_example_blast_output()
        display_result(blast_records)
        # TODO(roshan) write some tests for your functions
        # https://docs.python.org/3/library/unittest.html
        self.assertTrue(1 == 1)
