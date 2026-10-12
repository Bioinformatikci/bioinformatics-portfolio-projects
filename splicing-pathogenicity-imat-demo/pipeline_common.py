"""Shared, minimal functions for the public one-cohort demonstration."""
from __future__ import annotations

import warnings
from pathlib import Path

import cobra
import numpy as np
import pandas as pd
from cobamp.wrappers.cobra import COBRAModelObjectReader
from cobra.util.array import create_stoichiometric_matrix
from scipy.stats import fisher_exact, norm
from troppo.methods.reconstruction.imat import IMAT, IMATProperties


ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data"
RAW = DATA / "raw"
RESULTS = ROOT / "results"
MODEL_PATH = ROOT / "model" / "Human-GEM-1.12.0.xml"
SEED = 20261010
PERMUTATIONS = 1000
CONDITIONS = ("expression_only", "splicing_informed", "splicing_pathogenicity_informed")
ESSENTIAL_BOUNDS = {"MAR09034": (None, -0.01), "MAR09048": (None, -0.01), "MAR13082": (0.0001, None)}


def require(path: Path) -> Path:
    if not path.exists():
        raise FileNotFoundError(f"Required file is missing: {path}. Run the preceding script first.")
    return path


def normalise_gene_id(value: object) -> str:
    return str(value).strip().removeprefix("G_").split(".", 1)[0]


def omitnan_min(values) -> float:
    """COBAMP-compatible AND: minimum among available gene scores."""
    present = [float(x) for x in values if pd.notna(x)]
    return min(present) if present else np.nan


def omitnan_max(values) -> float:
    """COBAMP-compatible OR: maximum among available gene scores."""
    present = [float(x) for x in values if pd.notna(x)]
    return max(present) if present else np.nan


def load_model(*, require_gurobi: bool = True):
    require(MODEL_PATH)
    if require_gurobi:
        cobra.Configuration().solver = "gurobi"
    model = cobra.io.read_sbml_model(str(MODEL_PATH))
    for reaction_id, (lower, upper) in ESSENTIAL_BOUNDS.items():
        reaction = model.reactions.get_by_id(reaction_id)
        if lower is not None:
            reaction.lower_bound = lower
        if upper is not None:
            reaction.upper_bound = upper
    return model


def condition_dir(condition: str) -> Path:
    output = RESULTS / condition
    output.mkdir(parents=True, exist_ok=True)
    return output


def quartiles(values: np.ndarray) -> tuple[float, float]:
    values = np.asarray(values, dtype=float).ravel()
    median = float(np.median(values))
    return float(np.median(values[values < median])), float(np.median(values[values > median]))


def prepare_expression() -> tuple[pd.DataFrame, pd.DataFrame, tuple[float, float]]:
    expression = pd.read_csv(require(RAW / "expression.tsv"), sep="\t")
    groups = pd.read_csv(require(RAW / "metadata.tsv"), sep="\t")
    if set(groups.columns) != {"sample_id", "group"} or not set(groups.group).issubset({"ad", "control"}):
        raise ValueError("metadata.tsv must have sample_id and group (ad/control).")
    expression = expression.loc[expression.biotype.eq("protein_coding")].copy()
    expression["gene_id"] = expression.gene_id.map(normalise_gene_id)
    samples = groups.sample_id.astype(str).tolist()
    clean = expression.set_index("gene_id")[samples].apply(pd.to_numeric, errors="raise")
    shift = max(0.0, -float(np.nanmin(clean.to_numpy(dtype=float))))
    clean = clean + shift
    clean = clean.groupby(level=0).mean()
    low, high = quartiles(clean.mean(axis=1).to_numpy())
    return clean, groups, (low, high)


def map_to_model(expression: pd.DataFrame, model) -> pd.DataFrame:
    genes = list(dict.fromkeys(normalise_gene_id(gene.id) for gene in model.genes))
    return expression.reindex([gene for gene in genes if gene in expression.index])


def zero_ad_targets(before: pd.DataFrame, groups: pd.DataFrame, targets: pd.DataFrame, category: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    after = before.copy()
    ad_samples = set(groups.loc[groups.group.eq("ad"), "sample_id"].astype(str))
    audit = []
    for row in targets.itertuples(index=False):
        sample, gene = str(row.sample_id), normalise_gene_id(row.gene_id)
        if sample not in ad_samples or sample not in after.columns or gene not in after.index:
            raise ValueError(f"Invalid {category} target: {sample}:{gene}")
        value = float(after.at[gene, sample])
        after.at[gene, sample] = 0.0
        audit.append({"sample_id": sample, "gene_id": gene, "category": category,
                      "value_before": value, "value_after": 0.0,
                      "status": "newly_zeroed" if value != 0 else "already_zero"})
    return after, pd.DataFrame(audit)


def validate_zeroing(before: pd.DataFrame, after: pd.DataFrame, groups: pd.DataFrame, targets: pd.DataFrame) -> pd.DataFrame:
    if not before.index.equals(after.index) or not before.columns.equals(after.columns):
        raise RuntimeError("Expression axes changed during zeroing.")
    allowed = np.zeros(before.shape, dtype=bool)
    row_number = {gene: i for i, gene in enumerate(before.index)}
    column_number = {sample: i for i, sample in enumerate(before.columns)}
    for target in targets.itertuples(index=False):
        allowed[row_number[normalise_gene_id(target.gene_id)], column_number[str(target.sample_id)]] = True
    changed = before.to_numpy() != after.to_numpy()
    controls = [column_number[s] for s in groups.loc[groups.group.eq("control"), "sample_id"].astype(str)]
    summary = pd.DataFrame([{
        "allowed_target_cells": int(allowed.sum()), "changed_cells": int(changed.sum()),
        "unauthorised_changed_cells": int((changed & ~allowed).sum()),
        "control_cells_changed": int(changed[:, controls].sum()),
        "zero_to_nonzero_cells": int(((before.to_numpy() == 0) & (after.to_numpy() != 0)).sum()),
        "targets_not_zero_after": int((after.to_numpy()[allowed] != 0).sum()),
    }])
    summary["passed"] = (summary.unauthorised_changed_cells.eq(0) & summary.control_cells_changed.eq(0)
                          & summary.zero_to_nonzero_cells.eq(0) & summary.targets_not_zero_after.eq(0))
    if not bool(summary.at[0, "passed"]):
        raise RuntimeError("Zeroing quality control failed.")
    return summary


def score_reactions(expression: pd.DataFrame, model) -> pd.DataFrame:
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="Could not normalize this rule: .*")
        reader = COBRAModelObjectReader(model, gpr_gene_parse_function=normalise_gene_id)
    scores = {sample: reader.get_reaction_scores(expression[sample].dropna().to_dict(), and_fx=omitnan_min, or_fx=omitnan_max)
              for sample in expression.columns}
    return pd.DataFrame(scores).replace({None: np.nan})


def solve_imat_matrix(scores: pd.DataFrame, groups: pd.DataFrame, model, cutoffs: tuple[float, float]) -> tuple[pd.DataFrame, pd.DataFrame]:
    reaction_ids = [reaction.id for reaction in model.reactions]
    scores = scores.reindex(reaction_ids)
    stoichiometry = create_stoichiometric_matrix(model, array_type="dense", dtype=np.float64)
    lower = np.array([reaction.lower_bound for reaction in model.reactions])
    upper = np.array([reaction.upper_bound for reaction in model.reactions])
    matrix = pd.DataFrame(0, index=reaction_ids, columns=groups.sample_id.astype(str), dtype=np.int8)
    diagnostics = []
    low_cutoff, high_cutoff = cutoffs
    for number, sample in enumerate(matrix.columns, 1):
        values = scores[sample].to_numpy(dtype=float)
        properties = IMATProperties(exp_vector=values, exp_thresholds=(low_cutoff, high_cutoff), tolerance=1e-8, epsilon=1)
        algorithm = IMAT(stoichiometry, lower, upper, properties)
        kept = np.asarray(algorithm.run(), dtype=int)
        status = algorithm.sol.status()
        if status != "optimal":
            raise RuntimeError(f"{sample}: TROPPO iMAT status={status}")
        matrix.iloc[kept, matrix.columns.get_loc(sample)] = 1
        diagnostics.append({"sample_id": sample, "group": groups.set_index("sample_id").at[sample, "group"],
                            "active_reactions": int(matrix[sample].sum()), "imat_objective": float(algorithm.sol.objective_value()),
                            "high_reactions": int(np.isfinite(values).sum() and (values >= high_cutoff).sum()),
                            "low_reactions": int(np.isfinite(values).sum() and ((values >= 0) & (values < low_cutoff)).sum()),
                            "troppo_status": status, "result_type": "reaction_matrix_no_gem"})
        print(f"{sample}: TROPPO iMAT {number}/{len(matrix.columns)}", flush=True)
    matrix.index.name = "reaction_id"
    return matrix, pd.DataFrame(diagnostics)


def run_fisher(matrix: pd.DataFrame, groups: pd.DataFrame) -> pd.DataFrame:
    ad = groups.loc[groups.group.eq("ad"), "sample_id"].astype(str).tolist()
    control = groups.loc[groups.group.eq("control"), "sample_id"].astype(str).tolist()
    rows = []
    for reaction_id, values in matrix.iterrows():
        ad_present, control_present = int(values[ad].sum()), int(values[control].sum())
        odds, p_value = fisher_exact([[ad_present, len(ad) - ad_present], [control_present, len(control) - control_present]])
        rows.append({"reaction_id": reaction_id, "ad_present": ad_present, "ad_absent": len(ad) - ad_present,
                     "control_present": control_present, "control_absent": len(control) - control_present,
                     "odds_ratio": odds, "p_value": p_value})
    return pd.DataFrame(rows).sort_values("p_value", kind="stable")


def reporter_like(fisher: pd.DataFrame, subsystem: pd.DataFrame) -> pd.DataFrame:
    """Permutation reporter-like statistic; deliberately no hypergeometric test."""
    joined = subsystem.merge(fisher[["reaction_id", "p_value"]], on="reaction_id", how="inner")
    joined["z_score"] = norm.isf(joined.p_value.clip(lower=np.finfo(float).tiny, upper=1 - 1e-15))
    rng = np.random.default_rng(SEED)
    all_z = joined.z_score.to_numpy()
    rows = []
    for pathway, part in joined.groupby("subsystem", sort=True):
        observed = float(part.z_score.mean())
        size = len(part)
        null = np.array([rng.choice(all_z, size=size, replace=False).mean() for _ in range(PERMUTATIONS)])
        p_value = (1 + int((null >= observed).sum())) / (PERMUTATIONS + 1)
        rows.append({"subsystem": pathway, "reaction_count": size, "reporter_score": observed,
                     "reporter_p_value": p_value})
    result = pd.DataFrame(rows).sort_values("reporter_p_value", kind="stable").reset_index(drop=True)
    result["reporter_p_adjusted_bh"] = _bh(result.reporter_p_value.to_numpy())
    return result


def _bh(p_values: np.ndarray) -> np.ndarray:
    order = np.argsort(p_values)
    ranked = p_values[order] * len(p_values) / np.arange(1, len(p_values) + 1)
    adjusted = np.minimum.accumulate(ranked[::-1])[::-1]
    result = np.empty_like(adjusted)
    result[order] = np.minimum(adjusted, 1.0)
    return result
