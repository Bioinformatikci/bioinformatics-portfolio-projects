"""Layer 2: zero AD splicing targets, then rerun COBAMP + TROPPO iMAT."""
from __future__ import annotations

import pandas as pd

from pipeline_common import (RAW, RESULTS, condition_dir, load_model, require, score_reactions,
                             solve_imat_matrix, validate_zeroing, zero_ad_targets)


def main() -> None:
    baseline = RESULTS / "expression_only" / "processed"
    expression = pd.read_csv(require(baseline / "expression_mapped.tsv"), sep="\t", index_col=0)
    groups = pd.read_csv(require(baseline / "sample_groups.tsv"), sep="\t")
    cutoff_table = pd.read_csv(require(baseline / "expression_thresholds.tsv"), sep="\t")
    cutoffs = tuple(cutoff_table.value.astype(float))
    targets = pd.read_csv(require(RAW / "splicing_targets.tsv"), sep="\t").drop_duplicates(["sample_id", "gene_id"])
    after, audit = zero_ad_targets(expression, groups, targets, "splicing")
    qc = validate_zeroing(expression, after, groups, targets)
    output = condition_dir("splicing_informed")
    processed, matrix_dir = output / "processed", output / "matrix"
    processed.mkdir(exist_ok=True); matrix_dir.mkdir(exist_ok=True)
    after.to_csv(processed / "expression_mapped_post_splicing.tsv", sep="\t")
    groups.to_csv(processed / "sample_groups.tsv", sep="\t", index=False)
    targets.to_csv(processed / "splicing_targets.tsv", sep="\t", index=False)
    audit.to_csv(processed / "splicing_zeroing_audit.tsv", sep="\t", index=False)
    qc.to_csv(processed / "splicing_zeroing_qc.tsv", sep="\t", index=False)
    model = load_model()
    scores = score_reactions(after, model)
    scores.to_csv(processed / "reaction_scores.tsv", sep="\t")
    matrix, diagnostics = solve_imat_matrix(scores, groups, model, cutoffs, matrix_dir / "checkpoints",
                                            RESULTS / "expression_only" / "matrix" / "checkpoints")
    matrix.to_csv(matrix_dir / "reaction_activity_matrix.tsv", sep="\t")
    diagnostics.to_csv(matrix_dir / "sample_summary.tsv", sep="\t", index=False)
    print(f"Splicing-informed reaction matrix written to {matrix_dir}")


if __name__ == "__main__":
    main()
