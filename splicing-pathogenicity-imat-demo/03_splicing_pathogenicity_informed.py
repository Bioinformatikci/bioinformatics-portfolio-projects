"""Layer 3: add strict pathogenicity targets to the post-splicing expression."""
from __future__ import annotations

import pandas as pd

from pipeline_common import (RAW, RESULTS, condition_dir, load_model, require, score_reactions,
                             solve_imat_matrix, validate_zeroing, zero_ad_targets)


def main() -> None:
    splice_dir = RESULTS / "splicing_informed" / "processed"
    baseline_dir = RESULTS / "expression_only" / "processed"
    expression = pd.read_csv(require(splice_dir / "expression_mapped_post_splicing.tsv"), sep="\t", index_col=0)
    groups = pd.read_csv(require(baseline_dir / "sample_groups.tsv"), sep="\t")
    cutoffs = tuple(pd.read_csv(require(baseline_dir / "expression_thresholds.tsv"), sep="\t").value.astype(float))
    ranks = pd.read_csv(require(RAW / "pathogenicity_ranks.tsv"), sep="\t")
    ad_samples = groups.loc[groups.group.eq("ad"), "sample_id"].astype(str).tolist()
    rows = []
    for sample in ad_samples:
        if sample not in ranks.columns:
            raise ValueError(f"Pathogenicity rank column is missing for {sample}.")
        selected = ranks.loc[pd.to_numeric(ranks[sample], errors="coerce").gt(0.95), ["gene_id", "gene_symbol"]].copy()
        selected.insert(0, "sample_id", sample)
        selected["pathogenicity_rank"] = pd.to_numeric(ranks.loc[selected.index, sample])
        rows.append(selected)
    targets = pd.concat(rows, ignore_index=True).drop_duplicates(["sample_id", "gene_id"])
    after, audit = zero_ad_targets(expression, groups, targets, "pathogenicity_rank_gt_0_95")
    qc = validate_zeroing(expression, after, groups, targets)
    splice_targets = pd.read_csv(require(splice_dir / "splicing_targets.tsv"), sep="\t")
    membership = pd.concat([
        splice_targets.assign(target_source="splicing")[["sample_id", "gene_id", "target_source"]],
        targets.assign(target_source="pathogenicity")[["sample_id", "gene_id", "target_source"]],
    ]).groupby(["sample_id", "gene_id"], as_index=False).target_source.agg(";".join)
    output = condition_dir("splicing_pathogenicity_informed")
    processed, matrix_dir = output / "processed", output / "matrix"
    processed.mkdir(exist_ok=True); matrix_dir.mkdir(exist_ok=True)
    after.to_csv(processed / "expression_mapped_post_splicing_pathogenicity.tsv", sep="\t")
    targets.to_csv(processed / "pathogenicity_targets_rank_gt_0_95.tsv", sep="\t", index=False)
    membership.to_csv(processed / "combined_target_membership.tsv", sep="\t", index=False)
    audit.to_csv(processed / "pathogenicity_zeroing_audit.tsv", sep="\t", index=False)
    qc.to_csv(processed / "pathogenicity_zeroing_qc.tsv", sep="\t", index=False)
    model = load_model()
    scores = score_reactions(after, model)
    scores.to_csv(processed / "reaction_scores.tsv", sep="\t")
    matrix, diagnostics = solve_imat_matrix(scores, groups, model, cutoffs)
    matrix.to_csv(matrix_dir / "reaction_activity_matrix.tsv", sep="\t")
    diagnostics.to_csv(matrix_dir / "sample_summary.tsv", sep="\t", index=False)
    print(f"Splicing + pathogenicity-informed reaction matrix written to {matrix_dir}")


if __name__ == "__main__":
    main()
