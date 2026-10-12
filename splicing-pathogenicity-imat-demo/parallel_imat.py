"""Optional two-worker runner for layers 2 and 3; original iMAT scripts stay usable.

Each worker solves a disjoint set of 50 AD-labelled samples with the same
COBAMP scores, Human-GEM model, and TROPPO IMAT call as the serial scripts.
Controls are copied only after their expression is verified unchanged. A layer
matrix is written only after all worker checkpoints pass fingerprint and axis
validation. The expression-only process is never managed by this script.
"""
from __future__ import annotations

import argparse
import gc
import hashlib
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from pipeline_common import (RAW, RESULTS, ROOT, load_model, require, score_fingerprint,
                             score_reactions, solve_imat_matrix, validate_zeroing, zero_ad_targets)


LAYERS = {
    "splicing_informed": ("expression_only", "expression_mapped_post_splicing.tsv"),
    "splicing_pathogenicity_informed": ("splicing_informed", "expression_mapped_post_splicing_pathogenicity.tsv"),
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def input_expression(condition: str) -> tuple[pd.DataFrame, pd.DataFrame, tuple[float, float]]:
    baseline = RESULTS / "expression_only" / "processed"
    groups = pd.read_csv(require(baseline / "sample_groups.tsv"), sep="\t")
    if groups.group.value_counts().to_dict() != {"ad": 100, "control": 100} or groups.sample_id.duplicated().any():
        raise ValueError("Expected 100 unique AD and 100 unique control sample IDs")
    cutoff_table = pd.read_csv(require(baseline / "expression_thresholds.tsv"), sep="\t")
    values = dict(zip(cutoff_table.metric, cutoff_table.value))
    cutoffs = (float(values["low_cutoff"]), float(values["high_cutoff"]))
    previous = LAYERS[condition][0]
    path = (baseline / "expression_mapped.tsv" if previous == "expression_only"
            else RESULTS / previous / "processed" / "expression_mapped_post_splicing.tsv")
    expression = pd.read_csv(require(path), sep="\t", index_col=0)
    if expression.columns.astype(str).tolist() != groups.sample_id.astype(str).tolist():
        raise ValueError("Expression and metadata sample order differ")
    if expression.index.duplicated().any():
        raise ValueError("Duplicate model gene IDs")
    return expression, groups, cutoffs


def prepare_splicing(before: pd.DataFrame, groups: pd.DataFrame, processed: Path) -> pd.DataFrame:
    targets = pd.read_csv(require(RAW / "splicing_targets.tsv"), sep="\t").drop_duplicates(["sample_id", "gene_id"])
    after, audit = zero_ad_targets(before, groups, targets, "splicing")
    qc = validate_zeroing(before, after, groups, targets)
    groups.to_csv(processed / "sample_groups.tsv", sep="\t", index=False)
    targets.to_csv(processed / "splicing_targets.tsv", sep="\t", index=False)
    audit.to_csv(processed / "splicing_zeroing_audit.tsv", sep="\t", index=False)
    qc.to_csv(processed / "splicing_zeroing_qc.tsv", sep="\t", index=False)
    return after


def prepare_pathogenicity(before: pd.DataFrame, groups: pd.DataFrame, processed: Path) -> pd.DataFrame:
    ranks = pd.read_csv(require(RAW / "pathogenicity_ranks.tsv"), sep="\t")
    ad = groups.loc[groups.group.eq("ad"), "sample_id"].astype(str).tolist()
    rows = []
    for sample in ad:
        if sample not in ranks.columns:
            raise ValueError(f"Missing pathogenicity rank column: {sample}")
        selected = ranks.loc[pd.to_numeric(ranks[sample], errors="coerce").gt(0.95),
                             ["gene_id", "gene_symbol"]].copy()
        selected.insert(0, "sample_id", sample)
        selected["pathogenicity_rank"] = pd.to_numeric(ranks.loc[selected.index, sample])
        rows.append(selected)
    targets = pd.concat(rows, ignore_index=True).drop_duplicates(["sample_id", "gene_id"])
    after, audit = zero_ad_targets(before, groups, targets, "pathogenicity_rank_gt_0_95")
    qc = validate_zeroing(before, after, groups, targets)
    splice_targets = pd.read_csv(require(RESULTS / "splicing_informed" / "processed" / "splicing_targets.tsv"), sep="\t")
    membership = pd.concat([
        splice_targets.assign(target_source="splicing")[["sample_id", "gene_id", "target_source"]],
        targets.assign(target_source="pathogenicity")[["sample_id", "gene_id", "target_source"]],
    ]).groupby(["sample_id", "gene_id"], as_index=False).target_source.agg(";".join)
    targets.to_csv(processed / "pathogenicity_targets_rank_gt_0_95.tsv", sep="\t", index=False)
    membership.to_csv(processed / "combined_target_membership.tsv", sep="\t", index=False)
    audit.to_csv(processed / "pathogenicity_zeroing_audit.tsv", sep="\t", index=False)
    qc.to_csv(processed / "pathogenicity_zeroing_qc.tsv", sep="\t", index=False)
    return after


def prepare(condition: str) -> tuple[pd.DataFrame, tuple[float, float]]:
    previous, expression_name = LAYERS[condition]
    require(RESULTS / previous / "matrix" / "reaction_activity_matrix.tsv")
    before, groups, cutoffs = input_expression(condition)
    processed = RESULTS / condition / "processed"
    checkpoints = RESULTS / condition / "matrix" / "checkpoints"
    processed.mkdir(parents=True, exist_ok=True)
    after = (prepare_splicing(before, groups, processed) if condition == "splicing_informed"
             else prepare_pathogenicity(before, groups, processed))
    control = groups.loc[groups.group.eq("control"), "sample_id"].astype(str).tolist()
    if not np.array_equal(before[control].to_numpy(), after[control].to_numpy()):
        raise ValueError("Control expression changed")
    target = processed / expression_name
    temp = target.with_suffix(".tsv.tmp")
    after.to_csv(temp, sep="\t")
    if checkpoints.exists() and any(checkpoints.glob("*.json")) and target.exists() and sha256(temp) != sha256(target):
        temp.unlink()
        raise ValueError("Prepared expression changed while sample checkpoints exist")
    temp.replace(target)
    return groups, cutoffs


def ad_split(groups: pd.DataFrame, part: int) -> list[str]:
    ad = groups.loc[groups.group.eq("ad"), "sample_id"].astype(str).tolist()
    if len(ad) != 100 or part not in (0, 1):
        raise ValueError("The two-worker runner requires 100 AD samples and part 0 or 1")
    return ad[part * 50:(part + 1) * 50]


def worker(condition: str, part: int) -> None:
    _, expression_name = LAYERS[condition]
    _, groups, cutoffs = input_expression(condition)
    samples = ad_split(groups, part)
    expression = pd.read_csv(require(RESULTS / condition / "processed" / expression_name), sep="\t", index_col=0)
    model = load_model()
    scores = score_reactions(expression[samples], model)
    shard_groups = groups.set_index("sample_id").loc[samples].reset_index()
    previous = LAYERS[condition][0]
    matrix_dir = RESULTS / condition / "matrix"
    _, diagnostics = solve_imat_matrix(scores, shard_groups, model, cutoffs,
                                        matrix_dir / "checkpoints",
                                        RESULTS / previous / "matrix" / "checkpoints")
    shard_dir = matrix_dir / "shards"
    shard_dir.mkdir(exist_ok=True)
    diagnostics.to_csv(shard_dir / f"part_{part}_summary.tsv", sep="\t", index=False)
    print(f"Part {part}: {len(samples)} AD reaction vectors complete", flush=True)


def memory_gate() -> None:
    path = Path("/proc/meminfo")
    if path.exists():
        entry = next(line for line in path.read_text().splitlines() if line.startswith("MemAvailable:"))
        available_kib = int(entry.split()[1])
        if available_kib < 14 * 1024 * 1024:
            raise RuntimeError("Two iMAT workers need at least 14 GiB available RAM; run a serial layer or free memory")


def verify_and_merge(condition: str, groups: pd.DataFrame, cutoffs: tuple[float, float]) -> None:
    previous, expression_name = LAYERS[condition]
    model = load_model()
    reaction_ids = [reaction.id for reaction in model.reactions]
    previous_dir = RESULTS / previous / "matrix"
    prior = pd.read_csv(require(previous_dir / "reaction_activity_matrix.tsv"), sep="\t", index_col=0)
    prior_summary = pd.read_csv(require(previous_dir / "sample_summary.tsv"), sep="\t")
    samples = groups.sample_id.astype(str).tolist()
    if (prior.shape != (13070, 200) or prior.index.tolist() != reaction_ids
            or prior.columns.tolist() != samples or not prior.isin([0, 1]).all().all()):
        raise ValueError("Previous layer matrix axes differ from verified Human-GEM and metadata")
    ad = groups.loc[groups.group.eq("ad"), "sample_id"].astype(str).tolist()
    control = groups.loc[groups.group.eq("control"), "sample_id"].astype(str).tolist()
    expression = pd.read_csv(require(RESULTS / condition / "processed" / expression_name), sep="\t", index_col=0)
    scores = score_reactions(expression[ad], model).reindex(reaction_ids)
    current_dir = RESULTS / condition / "matrix"
    shard = pd.concat([pd.read_csv(require(current_dir / "shards" / f"part_{i}_summary.tsv"), sep="\t")
                       for i in (0, 1)], ignore_index=True)
    if shard.sample_id.astype(str).tolist() != ad or shard.sample_id.duplicated().any():
        raise ValueError("Worker sample lists overlap or are out of order")
    matrix = prior.copy()
    for sample in ad:
        payload = json.loads(require(current_dir / "checkpoints" / f"{sample}.json").read_text(encoding="utf-8"))
        expected = score_fingerprint(scores[sample].to_numpy(dtype=float), cutoffs)
        kept = np.asarray(payload["kept_indices"], dtype=int)
        if (payload.get("sample_id") != sample or payload.get("fingerprint") != expected
                or payload.get("reaction_count") != len(reaction_ids)
                or payload.get("solver", "GUROBI") != "GUROBI"
                or payload.get("troppo_status") != "optimal"
                or len(np.unique(kept)) != len(kept) or (kept < 0).any() or (kept >= len(reaction_ids)).any()):
            raise ValueError(f"Invalid or stale iMAT checkpoint: {sample}")
        matrix[sample] = 0
        matrix.iloc[kept, matrix.columns.get_loc(sample)] = 1
    if not matrix[control].equals(prior[control]):
        raise ValueError("Control reaction states differ from previous layer")
    for sample in control:
        original = require(previous_dir / "checkpoints" / f"{sample}.json")
        destination = current_dir / "checkpoints" / f"{sample}.json"
        shutil.copyfile(original, destination)
    old_control = prior_summary.set_index("sample_id").loc[control].reset_index()
    old_control["solver"] = "GUROBI"
    old_control["source"] = "reused_unchanged_control"
    summary = pd.concat([shard, old_control], ignore_index=True).set_index("sample_id").loc[samples].reset_index()
    if not np.array_equal(summary.active_reactions.to_numpy(), matrix.sum(axis=0).to_numpy()):
        raise ValueError("Merged sample summaries disagree with reaction matrix")
    if not summary.troppo_status.eq("optimal").all() or not summary.solver.eq("GUROBI").all():
        raise ValueError("Merged matrix contains non-optimal or non-Gurobi solutions")
    matrix.index.name = "reaction_id"
    matrix.to_csv(current_dir / "reaction_activity_matrix.tsv", sep="\t")
    summary.to_csv(current_dir / "sample_summary.tsv", sep="\t", index=False)
    print(f"Merged 100 AD vectors and 100 unchanged controls into {current_dir}")


def run_parent(condition: str) -> None:
    groups, cutoffs = prepare(condition)
    memory_gate()
    gc.collect()
    shard_dir = RESULTS / condition / "matrix" / "shards"
    shard_dir.mkdir(parents=True, exist_ok=True)
    processes = []
    log_handles = []
    try:
        for part in (0, 1):
            log = (shard_dir / f"part_{part}.log").open("w", encoding="utf-8")
            log_handles.append(log)
            command = [sys.executable, str(Path(__file__).resolve()), condition, "--worker-part", str(part)]
            processes.append(subprocess.Popen(command, cwd=ROOT, stdout=log, stderr=subprocess.STDOUT))
        print(f"Started disjoint 50+50 AD Gurobi workers for {condition}", flush=True)
        codes = [process.wait() for process in processes]
        if codes != [0, 0]:
            raise RuntimeError(f"A worker failed ({codes}); inspect logs under {shard_dir}. Checkpoints are retained")
    finally:
        for log in log_handles:
            log.close()
    verify_and_merge(condition, groups, cutoffs)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("condition", choices=LAYERS)
    parser.add_argument("--worker-part", type=int, choices=(0, 1), help=argparse.SUPPRESS)
    parser.add_argument("--plan", action="store_true", help="show the disjoint sample assignments without solving")
    parser.add_argument("--merge-only", action="store_true",
                        help="validate completed worker checkpoints and merge without solving")
    args = parser.parse_args()
    if args.plan:
        _, groups, _ = input_expression(args.condition)
        for part in (0, 1):
            samples = ad_split(groups, part)
            print(f"part {part}: {samples[0]}..{samples[-1]} ({len(samples)} AD samples)")
        return
    if os.environ.get("IMAT_SOLVER", "GUROBI").upper() != "GUROBI":
        raise ValueError("The published two-worker run requires IMAT_SOLVER=GUROBI")
    if args.merge_only:
        if args.worker_part is not None:
            parser.error("--merge-only cannot be combined with --worker-part")
        groups, cutoffs = input_expression(args.condition)[1:]
        verify_and_merge(args.condition, groups, cutoffs)
    elif args.worker_part is not None:
        worker(args.condition, args.worker_part)
    else:
        run_parent(args.condition)


if __name__ == "__main__":
    main()
