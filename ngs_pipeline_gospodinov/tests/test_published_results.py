"""Lightweight consistency checks for the published educational example outputs."""
import csv
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "variant_pipeline"

class PublishedResultsTests(unittest.TestCase):
    def test_variant_count_and_classes(self):
        with (OUT / "variants.tsv").open() as f:
            rows = list(csv.reader(f, delimiter="\t"))
        self.assertEqual(len(rows), 19)
        self.assertTrue(all(len(row) == 7 for row in rows))
        snps = sum(len(row[2]) == len(row[3]) == 1 for row in rows)
        self.assertEqual(snps, 12)
        self.assertEqual(len(rows)-snps, 7)

    def test_coverage_arithmetic(self):
        with (OUT / "coverage_summary.tsv").open() as f:
            rows = dict(csv.reader(f, delimiter="\t"))
        n = int(rows["Reference length (bp)"])
        self.assertEqual(n, 197394)
        self.assertEqual(int(rows["Bases covered >=1x"]) + int(rows["Bases below 1x"]), n)
        self.assertEqual(int(rows["Bases covered >=10x"]) + int(rows["Bases below 10x"]), n)
        self.assertAlmostEqual(100*int(rows["Bases covered >=10x"])/n, float(rows["Breadth >=10x (%)"]), places=3)

if __name__ == "__main__":
    unittest.main()
