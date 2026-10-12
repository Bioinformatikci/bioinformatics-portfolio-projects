"""Create a deterministic, non-biological one-cohort input example."""
from __future__ import annotations

import argparse
import shutil

import numpy as np
import pandas as pd

from pipeline_common import DATA, MODEL_PATH, RAW, load_model, normalise_gene_id


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--overwrite", action="store_true", help="replace generated synthetic data")
    args = parser.parse_args()
    if DATA.exists():
        if not args.overwrite:
            raise FileExistsError("data/ already exists; use --overwrite to regenerate the synthetic inputs.")
        shutil.rmtree(DATA)
    RAW.mkdir(parents=True)

    # Input creation needs the model gene/reaction IDs, but no MILP solver.
    model = load_model(require_gurobi=False)
    genes = list(dict.fromkeys(normalise_gene_id(gene.id) for gene in model.genes))
    if len(genes) < 20:
        raise RuntimeError("The Human-GEM model did not contain the expected gene set.")
    samples = ["DEMO_AD_01", "DEMO_AD_02", "DEMO_AD_03", "DEMO_AD_04",
               "DEMO_CONTROL_01", "DEMO_CONTROL_02", "DEMO_CONTROL_03", "DEMO_CONTROL_04"]
    groups = pd.DataFrame({"sample_id": samples, "group": ["ad"] * 4 + ["control"] * 4})
    groups.to_csv(RAW / "metadata.tsv", sep="\t", index=False)

    rng = np.random.default_rng(20261010)
    values = rng.gamma(shape=4.0, scale=2.0, size=(len(genes), len(samples)))
    expression = pd.DataFrame(values, columns=samples)
    expression.insert(0, "biotype", "protein_coding")
    expression.insert(0, "gene_id", genes)
    expression.to_csv(RAW / "expression.tsv", sep="\t", index=False)

    splice_rows = []
    pathogenicity_rows = []
    for sample_number, sample in enumerate(samples[:4]):
        splice_rows.extend([
            {"sample_id": sample, "gene_id": genes[10 + sample_number], "gene_name": f"Synthetic splice target {sample_number + 1}"},
            {"sample_id": sample, "gene_id": genes[20 + sample_number], "gene_name": f"Synthetic splice target {sample_number + 5}"},
        ])
        pathogenicity_rows.append((sample, genes[10 + sample_number], 0.99))  # overlaps splice
        pathogenicity_rows.append((sample, genes[30 + sample_number], 0.97))  # pathogenicity-only
    pd.DataFrame(splice_rows).to_csv(RAW / "splicing_targets.tsv", sep="\t", index=False)

    ranks = pd.DataFrame({"gene_id": genes, "gene_symbol": [f"GENE_{i + 1}" for i in range(len(genes))]})
    for sample in samples[:4]:
        ranks[sample] = 0.0
    for sample, gene, rank in pathogenicity_rows:
        ranks.loc[ranks.gene_id.eq(gene), sample] = rank
    ranks.to_csv(RAW / "pathogenicity_ranks.tsv", sep="\t", index=False)

    subsystem = pd.DataFrame({
        "reaction_id": [reaction.id for reaction in model.reactions],
        "subsystem": [reaction.subsystem or "Unannotated" for reaction in model.reactions],
    })
    subsystem.to_csv(RAW / "reaction_subsystem.tsv", sep="\t", index=False)
    pd.DataFrame([{
        "dataset": "synthetic_demo_cohort", "model": MODEL_PATH.name, "seed": 20261010,
        "ad_samples": 4, "control_samples": 4, "model_genes": len(genes), "model_reactions": len(model.reactions),
        "note": "Synthetic inputs only; not biological data.",
    }]).to_csv(DATA / "input_manifest.tsv", sep="\t", index=False)
    print(f"Created synthetic inputs for {len(samples)} samples and {len(genes)} Human-GEM genes in {RAW}")


if __name__ == "__main__":
    main()
