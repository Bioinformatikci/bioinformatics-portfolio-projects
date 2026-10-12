"""Lightweight checks for the committed synthetic inputs and model index."""
from __future__ import annotations

import csv
import gzip
import hashlib
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "synthetic_single_cohort"


class InputContract(unittest.TestCase):
    def test_supplied_input_hashes(self):
        with (SOURCE / "input_checksums.tsv").open(newline="", encoding="utf-8") as stream:
            expected = list(csv.DictReader(stream, delimiter="\t"))
        self.assertEqual({row["file"] for row in expected}, {
            "expression.tsv", "metadata.tsv", "pathogenicity_ranks.tsv",
            "splicing_targets.tsv", "synthetic_data_manifest.tsv"})
        for row in expected:
            name = row["file"]
            compressed = SOURCE / f"{name}.gz"
            path = compressed if compressed.exists() else SOURCE / name
            self.assertTrue(path.is_file(), name)
            digest = hashlib.sha256()
            opener = gzip.open if path.suffix == ".gz" else open
            with opener(path, "rb") as input_stream:
                for block in iter(lambda: input_stream.read(1024 * 1024), b""):
                    digest.update(block)
            self.assertEqual(digest.hexdigest(), row["sha256_uncompressed"], name)

    def test_sample_and_target_contract(self):
        with (SOURCE / "metadata.tsv").open(newline="", encoding="utf-8") as stream:
            metadata = list(csv.DictReader(stream, delimiter="\t"))
        self.assertEqual(len(metadata), 200)
        groups = {row["sample_id"]: row["group"] for row in metadata}
        self.assertEqual(len(groups), 200)
        self.assertEqual(sum(group == "ad" for group in groups.values()), 100)
        self.assertEqual(sum(group == "control" for group in groups.values()), 100)
        with (SOURCE / "splicing_targets.tsv").open(newline="", encoding="utf-8") as stream:
            targets = list(csv.DictReader(stream, delimiter="\t"))
        keys = [(row["sample_id"], row["gene_id"]) for row in targets]
        self.assertEqual(len(keys), len(set(keys)))
        self.assertEqual(len(keys), 800)
        self.assertTrue(all(groups[sample] == "ad" for sample, _ in keys))
        self.assertEqual(len({sample for sample, _ in keys}), 100)

    def test_verified_model_index_dimensions(self):
        path = ROOT / "reference" / "model_index.tsv.gz"
        self.assertTrue(path.is_file())
        with gzip.open(path, "rt", newline="", encoding="utf-8") as stream:
            rows = list(csv.DictReader(stream, delimiter="\t"))
        self.assertEqual(len(rows), 13070)
        self.assertEqual(len({row["reaction_id"] for row in rows}), 13070)
        self.assertEqual(sum(len(row["subsystems"].split(";")) for row in rows), 13071)


if __name__ == "__main__":
    unittest.main()
