"""Build a compact reaction/gene/subsystem index from Human-GEM v1.12.0.

The source SBML is verified before parsing. The index records associations,
not fluxes, inferred reaction activity, or iMAT results.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import xml.etree.ElementTree as ET
from pathlib import Path


ROOT = Path(__file__).resolve().parent
EXPECTED_SHA256 = "53a74c4495ca4af1334ce67b74cc09fc1ab30c7648eef561f478cac2363f13ab"
FBC = "{http://www.sbml.org/sbml/level3/version1/fbc/version2}"
GROUPS = "{http://www.sbml.org/sbml/level3/version1/groups/version1}"


def verify_model(path: Path) -> None:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(block)
    if digest.hexdigest() != EXPECTED_SHA256:
        raise ValueError(f"Human-GEM v1.12.0 SHA-256 mismatch: {path}")


def build_index(model: Path, output: Path) -> tuple[int, int]:
    verify_model(model)
    root = ET.parse(model).getroot()
    genes = {
        node.attrib[FBC + "id"]
        for node in root.iter()
        if node.tag.endswith("}geneProduct") and FBC + "id" in node.attrib
    }
    subsystem_by_reaction: dict[str, set[str]] = {}
    for group in root.iter(GROUPS + "group"):
        name = group.attrib.get(GROUPS + "name", "Unannotated")
        for member in group.iter(GROUPS + "member"):
            reaction = member.attrib.get(GROUPS + "idRef")
            if reaction:
                subsystem_by_reaction.setdefault(reaction, set()).add(name)

    rows = []
    for reaction in root.iter():
        if not reaction.tag.endswith("}reaction"):
            continue
        reaction_id = reaction.attrib["id"]
        linked_genes = {
            node.attrib[FBC + "geneProduct"]
            for node in reaction.iter()
            if node.tag.endswith("}geneProductRef") and FBC + "geneProduct" in node.attrib
        }
        if not linked_genes.issubset(genes):
            raise ValueError(f"Unknown model gene in {reaction_id}")
        rows.append((reaction_id, ";".join(sorted(subsystem_by_reaction.get(reaction_id, {"Unannotated"}))),
                     ";".join(sorted(linked_genes))))
    if len(genes) != 3067 or len(rows) != 13070:
        raise ValueError(f"Unexpected Human-GEM dimensions: {len(genes)} genes, {len(rows)} reactions")

    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("wb") as raw:
        with gzip.GzipFile(fileobj=raw, mode="wb", filename="", mtime=0, compresslevel=9) as compressed:
            with io.TextIOWrapper(compressed, encoding="utf-8", newline="") as text:
                writer = csv.writer(text, delimiter="\t", lineterminator="\n")
                writer.writerow(("reaction_id", "subsystems", "gene_ids"))
                writer.writerows(sorted(rows))
    return len(genes), len(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--model", type=Path, default=ROOT / "reference" / "Human-GEM_1.12.0.xml")
    parser.add_argument("--output", type=Path, default=ROOT / "reference" / "model_index.tsv.gz")
    args = parser.parse_args()
    genes, reactions = build_index(args.model, args.output)
    print(f"Indexed {genes} genes and {reactions} reactions from verified Human-GEM v1.12.0")


if __name__ == "__main__":
    main()
