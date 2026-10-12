"""Stage the supplied 100 AD / 100 control synthetic cohort for the iMAT scripts."""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import shutil
from pathlib import Path

from pipeline_common import DATA, RAW, ROOT


SOURCE = ROOT / "synthetic_single_cohort"
INDEX = ROOT / "reference" / "model_index.tsv.gz"
FILES = ("expression.tsv", "metadata.tsv", "splicing_targets.tsv", "pathogenicity_ranks.tsv")


def source_file(name: str) -> Path:
    compressed = SOURCE / f"{name}.gz"
    plain = SOURCE / name
    if compressed.exists():
        return compressed
    if plain.exists():
        return plain
    raise FileNotFoundError(f"Missing synthetic cohort input: {name}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--overwrite", action="store_true", help="replace an existing local data/ directory")
    args = parser.parse_args()
    if DATA.exists():
        if not args.overwrite:
            raise FileExistsError("data/ already exists; pass --overwrite only if it may be replaced.")
        shutil.rmtree(DATA)
    RAW.mkdir(parents=True)

    with (SOURCE / "input_checksums.tsv").open(newline="", encoding="utf-8") as stream:
        checksums = {row["file"]: row["sha256_uncompressed"] for row in csv.DictReader(stream, delimiter="\t")}
    for name in FILES:
        origin = source_file(name)
        digest = hashlib.sha256()
        opener = gzip.open if origin.suffix == ".gz" else open
        with opener(origin, "rb") as input_stream, (RAW / name).open("wb") as output_stream:
            for block in iter(lambda: input_stream.read(1024 * 1024), b""):
                digest.update(block)
                output_stream.write(block)
        if digest.hexdigest() != checksums[name]:
            raise ValueError(f"Checksum mismatch for {name}")

    if not INDEX.exists():
        raise FileNotFoundError("reference/model_index.tsv.gz is missing; run build_model_index.py.")
    with gzip.open(INDEX, "rt", newline="", encoding="utf-8") as stream, \
            (RAW / "reaction_subsystem.tsv").open("w", newline="", encoding="utf-8") as output:
        writer = csv.writer(output, delimiter="\t", lineterminator="\n")
        writer.writerow(("reaction_id", "SUBSYSTEM"))
        for row in csv.DictReader(stream, delimiter="\t"):
            for subsystem in row["subsystems"].split(";"):
                writer.writerow((row["reaction_id"], subsystem))
    shutil.copyfile(SOURCE / "synthetic_data_manifest.tsv", DATA / "input_manifest.tsv")
    print(f"Staged the supplied synthetic cohort in {RAW}")


if __name__ == "__main__":
    main()
