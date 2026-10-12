"""Checks for the committed, dependency-free synthetic showcase."""
from __future__ import annotations

import csv
import hashlib
import importlib.util
import re
import tempfile
import unittest
import xml.etree.ElementTree as ET
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SHOWCASE = ROOT / "showcase"
SPEC = importlib.util.spec_from_file_location("preview_demo", ROOT / "preview_demo.py")
assert SPEC and SPEC.loader
preview = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(preview)


def table(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


class ShowcaseTests(unittest.TestCase):
    def test_committed_files_are_reproducible(self) -> None:
        with tempfile.TemporaryDirectory(prefix="preview-check-", dir=ROOT) as directory:
            generated = Path(directory)
            preview.generate(generated)
            committed_paths = sorted(path.relative_to(SHOWCASE) for path in SHOWCASE.rglob("*") if path.is_file())
            generated_paths = sorted(path.relative_to(generated) for path in generated.rglob("*") if path.is_file())
            self.assertEqual(committed_paths, generated_paths)
            for relative in committed_paths:
                self.assertEqual(hashlib.sha256((SHOWCASE / relative).read_bytes()).digest(),
                                 hashlib.sha256((generated / relative).read_bytes()).digest(), str(relative))

    def test_only_allowed_ad_cells_change_and_rank_boundary_is_strict(self) -> None:
        raw = SHOWCASE / "data" / "raw"
        metadata = table(raw / "metadata.tsv")
        controls = {row["sample_id"] for row in metadata if row["group"] == "control"}
        ad = {row["sample_id"] for row in metadata if row["group"] == "ad"}
        self.assertEqual((len(ad), len(controls)), (4, 4))
        audit = table(SHOWCASE / "audits" / "zeroing.tsv")
        self.assertTrue(audit)
        self.assertTrue(all(row["sample_id"] in ad and row["value_after"] == "0.00" for row in audit))
        self.assertFalse(any(row["gene_id"] == "SYN_G010" for row in audit))
        self.assertEqual(sum(row["status"] == "already_zero" for row in audit), 4)
        matrices = [table(SHOWCASE / "results" / condition / "reaction_activity_matrix.tsv")
                    for condition in preview.CONDITIONS]
        for control in controls:
            for before, after in zip(matrices[0], matrices[-1]):
                self.assertEqual(before[control], after[control])

    def test_report_links_and_figures_are_valid(self) -> None:
        html = (SHOWCASE / "index.html").read_text(encoding="utf-8")
        self.assertIn("Synthetic preview only", html)
        for link in re.findall(r'(?:href|src)="([^"]+)"', html):
            self.assertTrue((SHOWCASE / link).is_file(), link)
        for figure in (SHOWCASE / "figures").glob("*.svg"):
            root = ET.parse(figure).getroot()
            self.assertEqual(root.tag, "{http://www.w3.org/2000/svg}svg")
            self.assertIsNotNone(root.find("{http://www.w3.org/2000/svg}title"))


if __name__ == "__main__":
    unittest.main()
