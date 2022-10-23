from pathlib import Path
from unittest import mock
from Bio.Blast import NCBIXML
import unittest
from flams.run_blast import _filter_blast
from flams.run_blast import ModificationHeader


TEST_BLAST_OUTPUT = Path(__file__).parent / "testfiles/blast_output.xml"


class FilterBlastTestCase(unittest.TestCase):
    def setUp(self) -> None:
        self.handle = open(TEST_BLAST_OUTPUT)
        self.blast_record = list(NCBIXML.parse(self.handle))[0]

    def tearDown(self) -> None:
        self.handle.close()

    def test_parse_title(self):
        title = "PLMD-7244|P25665|304 Acetylation [Escherichia coli (strain K12)]"
        mod = ModificationHeader.parse(title)
        self.assertEqual(mod.plmd_id, "PLMD-7244")
        self.assertEqual(mod.uniprot_id, "P25665")
        self.assertEqual(mod.position, 304)
        self.assertEqual(mod.modification, "Acetylation")
        self.assertEqual(mod.species, "Escherichia coli (strain K12)")

    def _count_hsps(self, blast_record, condition=lambda e: True):
        return len(
            [hsp for a in blast_record.alignments for hsp in a.hsps if condition(hsp)]
        )

    def test_filters_by_evalue(self):
        # given
        evalue = 1e-18

        all_hsps = self._count_hsps(self.blast_record)
        self.assertEqual(all_hsps, 112)

        expected_count = self._count_hsps(
            self.blast_record, lambda e: e.expect < evalue
        )
        self.assertEqual(expected_count, 98)

        # when
        # patch function so evalue is the only filter
        with mock.patch("flams.run_blast.does_hsp_match_query", return_value=True) as _:
            filtered = _filter_blast(self.blast_record, 1, 1, evalue)

        # then
        self.assertEqual(self._count_hsps(filtered), expected_count)
