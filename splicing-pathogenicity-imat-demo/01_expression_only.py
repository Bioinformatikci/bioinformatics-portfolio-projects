"""Layer 1: expression-only COBAMP scoring and TROPPO iMAT reaction matrix."""
from __future__ import annotations

import pandas as pd

from pipeline_common import (RESULTS, condition_dir, load_model, map_to_model, prepare_expression,
                             score_reactions, solve_imat_matrix)


def main() -> None:
    model = load_model()
    expression, groups, cutoffs = prepare_expression()
    mapped = map_to_model(expression, model)
    output = condition_dir("expression_only")
    processed, matrix_dir = output / "processed", output / "matrix"
    processed.mkdir(exist_ok=True); matrix_dir.mkdir(exist_ok=True)
    expression.to_csv(processed / "expression_clean.tsv", sep="\t")
    mapped.to_csv(processed / "expression_mapped.tsv", sep="\t")
    groups.to_csv(processed / "sample_groups.tsv", sep="\t", index=False)
    pd.DataFrame({"metric": ["low_cutoff", "high_cutoff"], "value": cutoffs}).to_csv(processed / "expression_thresholds.tsv", sep="\t", index=False)
    scores = score_reactions(mapped, model)
    scores.to_csv(processed / "reaction_scores.tsv", sep="\t")
    labels = scores.applymap(lambda x: "no_expression_evidence" if pd.isna(x) else "low" if x <= cutoffs[0] else "high" if x >= cutoffs[1] else "neutral")
    labels.to_csv(processed / "reaction_imat_labels.tsv", sep="\t")
    matrix, diagnostics = solve_imat_matrix(scores, groups, model, cutoffs, matrix_dir / "checkpoints")
    matrix.to_csv(matrix_dir / "reaction_activity_matrix.tsv", sep="\t")
    diagnostics.to_csv(matrix_dir / "sample_summary.tsv", sep="\t", index=False)
    print(f"Expression-only reaction matrix written to {matrix_dir}")


if __name__ == "__main__":
    main()
