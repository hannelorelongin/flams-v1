from pathlib import Path
import unittest
import csv

from Bio.Blast import NCBIXML
from flams.display import display_result


EXAMPLE_BLAST_OUTPUT = Path(__file__).parent / "testfiles/blast_output_P57703.xml"


class DisplayTestCase(unittest.TestCase):
    def test_display_simple(self):
        with EXAMPLE_BLAST_OUTPUT.open() as f:
            blast_records = list(NCBIXML.parse(f))
            display_result(
                output_filename="out_example.tsv", blast_records=blast_records
            )
        self.assertTrue(Path("out_example.tsv").is_file())
        out_example = Path("out_example.tsv")
        with out_example.open() as t:
            example_table = csv.reader(t, delimiter="\t")
            row1 = next(example_table)
            row2 = next(example_table)
            self.assertEqual(
                row2[0],
                "P25665",
                "The Uniprot Id is not correct in the example output file",
            )
            self.assertEqual(
                row2[1],
                "Escherichia coli (strain K12)",
                "The name of the species is not displayed correctly",
            )
            self.assertEqual(
                row2[2],
                "Acetylation",
                "The lysine modification is not displayed correctly",
            )
            self.assertEqual(
                int(row2[3]), 87, "The lysine location is not correctly displayed"
            )

            self.assertEqual(
                row2[4],
                "82-ARHQNKDGSVD-93",
                "The lysine window is not correctly displayed",
            )


if __name__ == "__main__":
    unittest.main()
