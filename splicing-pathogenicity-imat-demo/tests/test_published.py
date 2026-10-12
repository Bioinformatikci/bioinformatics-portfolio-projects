"""Check that every published figure and table is backed by complete matrices."""
from __future__ import annotations

import csv
import gzip
import hashlib
import json
import unittest
from html.parser import HTMLParser
from pathlib import Path
from urllib.parse import unquote, urlsplit


ROOT = Path(__file__).resolve().parents[1]
SHOWCASE = ROOT / "showcase"
CONDITIONS = ("expression_only", "splicing_informed", "splicing_pathogenicity_informed")
FIGURES = ("active_reactions.svg", "paired_model_size.svg",
           "target_expression.svg", "reporter_pathways.svg")


class LocalReferences(HTMLParser):
    def __init__(self):
        super().__init__()
        self.references = []

    def handle_starttag(self, tag, attrs):
        for name, value in attrs:
            if name in {"href", "src"} and value:
                self.references.append(value)


class PublishedContract(unittest.TestCase):
    def test_manifest_and_local_links(self):
        self.assertTrue((SHOWCASE / "index.html").is_file())
        with (SHOWCASE / "manifest.sha256").open(encoding="utf-8") as stream:
            rows = [line.rstrip("\n").split("  ", 1) for line in stream if line.strip()]
        listed = {name for _, name in rows}
        actual = {path.relative_to(SHOWCASE).as_posix() for path in SHOWCASE.rglob("*")
                  if path.is_file() and path.name != "manifest.sha256"}
        self.assertEqual(listed, actual)
        for expected, name in rows:
            self.assertEqual(hashlib.sha256((SHOWCASE / name).read_bytes()).hexdigest(), expected, name)
        parser = LocalReferences()
        parser.feed((SHOWCASE / "index.html").read_text(encoding="utf-8"))
        for reference in parser.references:
            parsed = urlsplit(reference)
            if parsed.scheme or parsed.netloc or not parsed.path:
                continue
            self.assertTrue((SHOWCASE / unquote(parsed.path)).is_file(), reference)
        for name in FIGURES:
            self.assertIn(f"figures/{name}", listed)

    def test_matrix_axes_controls_and_diagnostics(self):
        with (ROOT / "synthetic_single_cohort" / "metadata.tsv").open(newline="", encoding="utf-8") as stream:
            metadata = list(csv.DictReader(stream, delimiter="\t"))
        samples = [row["sample_id"] for row in metadata]
        control_indexes = [i for i, row in enumerate(metadata) if row["group"] == "control"]
        baseline_reactions = baseline_controls = None
        for condition in CONDITIONS:
            folder = SHOWCASE / "results" / condition
            with gzip.open(folder / "reaction_activity_matrix.tsv.gz", "rt", newline="", encoding="utf-8") as stream:
                reader = csv.reader(stream, delimiter="\t")
                header = next(reader)
                self.assertEqual(header[1:], samples, condition)
                sums = [0] * 200
                reactions, controls = [], []
                for row in reader:
                    self.assertEqual(len(row), 201, condition)
                    self.assertTrue(all(value in {"0", "1"} for value in row[1:]), condition)
                    reactions.append(row[0])
                    controls.append("".join(row[i + 1] for i in control_indexes))
                    for i, value in enumerate(row[1:]):
                        sums[i] += value == "1"
            self.assertEqual(len(reactions), 13070, condition)
            if baseline_reactions is None:
                baseline_reactions, baseline_controls = reactions, controls
            else:
                self.assertEqual(reactions, baseline_reactions, condition)
                self.assertEqual(controls, baseline_controls, condition)
            with (folder / "sample_summary.tsv").open(newline="", encoding="utf-8") as stream:
                summary = list(csv.DictReader(stream, delimiter="\t"))
            self.assertEqual([row["sample_id"] for row in summary], samples, condition)
            self.assertEqual([int(row["active_reactions"]) for row in summary], sums, condition)
            self.assertTrue(all(row["troppo_status"] == "optimal" and row["solver"] == "GUROBI"
                                for row in summary), condition)

    def test_statistics_and_metrics_contract(self):
        metrics = json.loads((SHOWCASE / "results" / "figure_metrics.json").read_text(encoding="utf-8"))
        self.assertIs(metrics["synthetic"], True)
        self.assertEqual(metrics["samples"], 200)
        target = metrics["splicing_targets"]
        self.assertEqual(target["target_pairs"], 800)
        self.assertEqual(target["low_pairs"] + target["neutral_pairs"] + target["high_pairs"], 800)
        self.assertEqual(metrics["paired_model_size_and_similarity"]["paired_ad_samples"], 100)
        for condition in CONDITIONS:
            folder = SHOWCASE / "results" / condition
            with gzip.open(folder / "fisher_all_reactions.tsv.gz", "rt", newline="", encoding="utf-8") as stream:
                fisher = list(csv.DictReader(stream, delimiter="\t"))
            self.assertEqual(len(fisher), 13070, condition)
            self.assertTrue({"reaction_id", "p_value", "p_adjusted_bh"}.issubset(fisher[0]), condition)
            with (folder / "reporter_like_all_pathways.tsv").open(newline="", encoding="utf-8") as stream:
                self.assertIn("reporter_p_corrected", next(csv.reader(stream, delimiter="\t")), condition)


if __name__ == "__main__":
    unittest.main()
