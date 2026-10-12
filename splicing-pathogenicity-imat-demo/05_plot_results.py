"""Create figures and a concise metrics report from the completed synthetic run."""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from pipeline_common import CONDITIONS, RAW, RESULTS, require


LABELS = {
    "expression_only": "Baseline expression",
    "splicing_informed": "Splicing adjusted",
    "splicing_pathogenicity_informed": "Splicing + rank adjusted",
}
NAVY, TEAL, ORANGE, MUTED = "#173349", "#2f7181", "#d6784c", "#60727e"


def style():
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11,
                         "svg.fonttype": "none", "axes.titleweight": "bold"})


def finish(fig, axes, path: Path, footnote: str, *, top: float = 0.94):
    for ax in np.atleast_1d(axes).ravel():
        if ax.axison:
            ax.spines[["top", "right"]].set_visible(False)
            ax.set_axisbelow(True)
            ax.grid(axis="y", color="#e3eaec", linewidth=0.8)
    fig.tight_layout(rect=(0.02, 0.11, 0.98, top))
    fig.text(0.04, 0.035, footnote, color=MUTED, fontsize=8.5)
    fig.savefig(path, bbox_inches="tight", facecolor="white")
    svg = path.read_text(encoding="utf-8")
    path.write_text("\n".join(line.rstrip() for line in svg.splitlines()) + "\n", encoding="utf-8")
    plt.close(fig)


def load_results() -> tuple[pd.DataFrame, dict]:
    groups = pd.read_csv(require(RAW / "metadata.tsv"), sep="\t")
    if groups.sample_id.duplicated().any() or groups.group.value_counts().to_dict() != {"ad": 100, "control": 100}:
        raise ValueError("Expected 100 unique AD-labelled and 100 unique control-labelled samples")
    samples = groups.sample_id.astype(str).tolist()
    result = {}
    reaction_ids = None
    for condition in CONDITIONS:
        root = RESULTS / condition
        matrix = pd.read_csv(require(root / "matrix" / "reaction_activity_matrix.tsv"), sep="\t", index_col=0)
        summary = pd.read_csv(require(root / "matrix" / "sample_summary.tsv"), sep="\t")
        fisher = pd.read_csv(require(root / "fisher" / "fisher_all_reactions.tsv"), sep="\t")
        reporter = pd.read_csv(require(root / "reporter_like" / "reporter_like_all_pathways.tsv"), sep="\t")
        if matrix.shape != (13070, 200) or list(matrix.columns) != samples or not matrix.isin([0, 1]).all().all():
            raise ValueError(f"Invalid reaction matrix: {condition}")
        if reaction_ids is None:
            reaction_ids = matrix.index.tolist()
        elif matrix.index.tolist() != reaction_ids:
            raise ValueError(f"Reaction order changed: {condition}")
        if summary.sample_id.astype(str).tolist() != samples or not np.array_equal(
                summary.active_reactions.to_numpy(), matrix.sum(axis=0).to_numpy()):
            raise ValueError(f"Sample summary does not match matrix: {condition}")
        if len(fisher) != 13070 or fisher.reaction_id.duplicated().any():
            raise ValueError(f"Incomplete Fisher table: {condition}")
        result[condition] = {"matrix": matrix, "summary": summary, "fisher": fisher, "reporter": reporter}
    control = groups.loc[groups.group.eq("control"), "sample_id"].astype(str).tolist()
    first = result[CONDITIONS[0]]["matrix"][control]
    for condition in CONDITIONS[1:]:
        if not first.equals(result[condition]["matrix"][control]):
            raise ValueError(f"Control reaction states changed in {condition}")
    return groups, result


def target_expression(groups: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    source = pd.read_csv(require(RAW / "splicing_targets.tsv"), sep="\t")
    audit = pd.read_csv(require(RESULTS / "splicing_informed" / "processed" / "splicing_zeroing_audit.tsv"), sep="\t")
    baseline = pd.read_csv(require(RESULTS / "expression_only" / "processed" / "expression_mapped.tsv"), sep="\t", index_col=0)
    cutoffs = dict(pd.read_csv(require(RESULTS / "expression_only" / "processed" / "expression_thresholds.tsv"), sep="\t").values)
    low, high = float(cutoffs["low_cutoff"]), float(cutoffs["high_cutoff"])
    if not low < high or source.duplicated(["sample_id", "gene_id"]).any():
        raise ValueError("Invalid splicing target IDs or expression cutoffs")
    ad = set(groups.loc[groups.group.eq("ad"), "sample_id"].astype(str))
    if not set(source.sample_id.astype(str)).issubset(ad):
        raise ValueError("Splicing targets include a non-AD sample")
    if set(zip(source.sample_id, source.gene_id)) != set(zip(audit.sample_id, audit.gene_id)):
        raise ValueError("Target list and zeroing audit differ")
    status = audit.set_index(["sample_id", "gene_id"])["status"]
    rows = []
    for row in source.itertuples(index=False):
        sample, gene = str(row.sample_id), str(row.gene_id)
        if gene not in baseline.index or sample not in baseline.columns:
            raise ValueError(f"Unmapped splicing target: {sample}/{gene}")
        value = float(baseline.at[gene, sample])
        category = "low" if value < low else "high" if value >= high else "neutral"
        rows.append({"sample_id": sample, "gene_id": gene, "expression_before": value,
                     "category": category, "zeroing_status": status.at[(sample, gene)]})
    frame = pd.DataFrame(rows).sort_values(["sample_id", "gene_id"], kind="stable")
    if frame.sample_id.nunique() != len(ad):
        raise ValueError("Some AD-labelled samples lack a splicing target")
    return frame, {"low_cutoff": low, "high_cutoff": high,
                   "target_pairs": len(frame), "samples_with_targets": frame.sample_id.nunique(),
                   "newly_zeroed_pairs": int(frame.zeroing_status.eq("newly_zeroed").sum()),
                   "already_zero_pairs": int(frame.zeroing_status.eq("already_zero").sum()),
                   "low_pairs": int(frame.category.eq("low").sum()),
                   "neutral_pairs": int(frame.category.eq("neutral").sum()),
                   "high_pairs": int(frame.category.eq("high").sum()),
                   "detected_splice_event_gene_count_available": False}


def overlap(left: set, right: set) -> dict:
    intersection = len(left & right)
    union = len(left | right)
    minimum = min(len(left), len(right))
    return {"left": len(left), "right": len(right), "intersection": intersection,
            "jaccard": intersection / union if union else None,
            "overlap_index": intersection / minimum if minimum else None}


def paired_metrics(groups: pd.DataFrame, result: dict) -> tuple[pd.DataFrame, dict]:
    ad = groups.loc[groups.group.eq("ad"), "sample_id"].astype(str).tolist()
    baseline = result[CONDITIONS[0]]["matrix"]
    splice = result[CONDITIONS[1]]["matrix"]
    final = result[CONDITIONS[2]]["matrix"]
    rows = []
    for sample in ad:
        before = set(baseline.index[baseline[sample].eq(1)])
        after = set(splice.index[splice[sample].eq(1)])
        similarity = overlap(before, after)
        rows.append({"sample_id": sample,
                     "baseline_active": len(before), "splicing_active": len(after),
                     "final_active": int(final[sample].sum()),
                     "baseline_minus_splicing": len(before) - len(after),
                     "splicing_minus_final": len(after) - int(final[sample].sum()),
                     "reaction_jaccard": similarity["jaccard"],
                     "reaction_overlap_index": similarity["overlap_index"]})
    frame = pd.DataFrame(rows)
    first_paths = set(result[CONDITIONS[0]]["reporter"].loc[
        result[CONDITIONS[0]]["reporter"].reporter_p_corrected.lt(0.05), "subsystem"])
    splice_paths = set(result[CONDITIONS[1]]["reporter"].loc[
        result[CONDITIONS[1]]["reporter"].reporter_p_corrected.lt(0.05), "subsystem"])
    series = frame.baseline_minus_splicing
    q1, q3 = (float(value) for value in series.quantile([0.25, 0.75]))
    result_metrics = {"paired_ad_samples": len(frame),
                      "baseline_minus_splicing_median": float(series.median()),
                      "baseline_minus_splicing_q1": q1,
                      "baseline_minus_splicing_q3": q3,
                      "baseline_minus_splicing_outliers": int(((series < q1 - 1.5 * (q3 - q1)) |
                                                              (series > q3 + 1.5 * (q3 - q1))).sum()),
                      "reaction_jaccard_median": float(frame.reaction_jaccard.median()),
                      "reaction_overlap_index_median": float(frame.reaction_overlap_index.median()),
                      "reporter_significant_pathway_overlap": overlap(first_paths, splice_paths)}
    return frame, result_metrics


def fisher_subsystem_table(result: dict) -> pd.DataFrame:
    annotation = pd.read_csv(require(RAW / "reaction_subsystem.tsv"), sep="\t")
    if annotation.columns.tolist() != ["reaction_id", "SUBSYSTEM"] or len(annotation) != 13071:
        raise ValueError("Unexpected Human-GEM subsystem annotation")
    rows = []
    for condition in CONDITIONS:
        fisher = result[condition]["fisher"]
        mapped = fisher.loc[fisher.p_value.lt(0.05), ["reaction_id"]].merge(annotation, on="reaction_id", how="left")
        significant = set(result[condition]["reporter"].loc[
            result[condition]["reporter"].reporter_p_corrected.lt(0.05), "subsystem"])
        for subsystem, part in mapped.groupby("SUBSYSTEM", dropna=False):
            rows.append({"condition": condition, "subsystem": subsystem,
                         "associated_fisher_reactions_raw_p_lt_0_05": part.reaction_id.nunique(),
                         "reporter_corrected_p_lt_0_05": subsystem in significant})
    if not rows:
        return pd.DataFrame(columns=["condition", "subsystem", "associated_fisher_reactions_raw_p_lt_0_05",
                                     "reporter_corrected_p_lt_0_05"])
    return pd.DataFrame(rows).sort_values(["condition", "associated_fisher_reactions_raw_p_lt_0_05"],
                                          ascending=[True, False], kind="stable")


def plot_activity(groups: pd.DataFrame, result: dict, out: Path):
    control = groups.loc[groups.group.eq("control"), "sample_id"].astype(str).tolist()
    ad = groups.loc[groups.group.eq("ad"), "sample_id"].astype(str).tolist()
    categories = [("Controls", CONDITIONS[0], control, TEAL),
                  ("Baseline AD", CONDITIONS[0], ad, ORANGE),
                  ("Splicing AD", CONDITIONS[1], ad, "#b85f45"),
                  ("Final AD", CONDITIONS[2], ad, "#865769")]
    fig, ax = plt.subplots(figsize=(10.4, 5.5))
    rng = np.random.default_rng(20261010)
    for i, (label, condition, samples, color) in enumerate(categories):
        counts = result[condition]["summary"].set_index("sample_id").loc[samples, "active_reactions"].to_numpy()
        ax.boxplot(counts, positions=[i], widths=0.42, patch_artist=True, showfliers=False,
                   medianprops={"color": NAVY, "linewidth": 2},
                   boxprops={"facecolor": color, "alpha": 0.25, "edgecolor": color},
                   whiskerprops={"color": color}, capprops={"color": color})
        ax.scatter(i + rng.normal(0, 0.065, len(counts)), counts, s=13, color=color, alpha=0.43, linewidths=0)
    ax.set_xticks(range(4), [item[0] for item in categories])
    ax.set_ylabel("Active reactions per sample")
    ax.set_title("Reaction activity across staged iMAT inputs", loc="left", color=NAVY, fontsize=16, pad=14)
    finish(fig, ax, out / "active_reactions.svg", "Synthetic cohort · n = 100 per group/stage · controls shown once because their input is unchanged")


def plot_paired_size(frame: pd.DataFrame, metrics: dict, out: Path):
    fig, ax = plt.subplots(figsize=(10.0, 5.3))
    series = [frame.baseline_minus_splicing.to_numpy(), frame.splicing_minus_final.to_numpy()]
    labels = ["Baseline − splicing", "Splicing − final"]
    colors = [TEAL, ORANGE]
    rng = np.random.default_rng(20261011)
    for i, (values, color) in enumerate(zip(series, colors)):
        ax.boxplot(values, positions=[i], widths=0.42, patch_artist=True, showfliers=False,
                   medianprops={"color": NAVY, "linewidth": 2},
                   boxprops={"facecolor": color, "alpha": 0.22, "edgecolor": color},
                   whiskerprops={"color": color}, capprops={"color": color})
        ax.scatter(i + rng.normal(0, 0.065, len(values)), values, s=13, color=color, alpha=0.48, linewidths=0)
    ax.axhline(0, color=NAVY, linewidth=1, linestyle="--", alpha=0.65)
    ax.set_xticks(range(2), labels)
    ax.set_ylabel("Before − after active reactions / paired AD sample")
    ax.set_title("Paired model-size change", loc="left", color=NAVY, fontsize=16, pad=14)
    finish(fig, ax, out / "paired_model_size.svg", "Positive values mean fewer active reactions after the indicated adjustment · all 100 AD pairs shown")


def plot_target_expression(frame: pd.DataFrame, metrics: dict, out: Path):
    fig, (left, right) = plt.subplots(1, 2, figsize=(10.8, 4.6), gridspec_kw={"width_ratios": [1.2, 1]})
    categories = [("low", "#7fa6b4"), ("neutral", "#d9a66b"), ("high", "#bd674e")]
    start = 0
    for category, color in categories:
        count = int(frame.category.eq(category).sum())
        left.barh([0], [count], left=start, color=color, height=0.46)
        if count:
            left.text(start + count / 2, 0, f"{category.title()}\n{count}",
                      ha="center", va="center", color=NAVY, weight="bold", fontsize=9)
        start += count
    left.set_yticks([0], ["Target pairs"])
    left.set_xlim(0, len(frame))
    left.set_xlabel("Gene × sample pairs")
    high_per_sample = frame.assign(is_high=frame.category.eq("high")).groupby("sample_id").is_high.sum()
    counts = high_per_sample.value_counts().sort_index()
    right.bar(counts.index, counts.values, width=0.75, color=TEAL)
    right.set_xlabel("High-expression targets / AD sample")
    right.set_ylabel("Number of AD samples")
    right.set_xticks(range(int(high_per_sample.max()) + 1))
    fig.suptitle("Target expression before splicing zeroing", x=0.04, ha="left", color=NAVY, fontsize=16, weight="bold")
    finish(fig, (left, right), out / "target_expression.svg",
           f"Descriptive gene-level categories using iMAT cutoffs: low < {metrics['low_cutoff']:.3f}; high ≥ {metrics['high_cutoff']:.3f} · {len(frame)} provided targets")


def plot_reporter(result: dict, out: Path):
    baseline = result[CONDITIONS[0]]["reporter"].set_index("subsystem")
    splice = result[CONDITIONS[1]]["reporter"].set_index("subsystem")
    significant = set(baseline.index[baseline.reporter_p_corrected.lt(0.05)]) | set(
        splice.index[splice.reporter_p_corrected.lt(0.05)])
    if not significant:
        fig, ax = plt.subplots(figsize=(9.8, 3.8))
        ax.axis("off")
        ax.text(0.5, 0.58, "No corrected reporter-like pathways pass p < 0.05",
                transform=ax.transAxes, ha="center", va="center", fontsize=15, color=NAVY)
        ax.text(0.5, 0.39, "The figure reflects the calculated synthetic outputs; no pathway was added for display.",
                transform=ax.transAxes, ha="center", color=MUTED)
        finish(fig, ax, out / "reporter_pathways.svg", "Archived R reporter-like method · 1,000 permutations · synthetic inputs")
        return
    ranked = sorted(significant, key=lambda name: min(
        baseline.at[name, "reporter_p_corrected"] if name in baseline.index else 1.0,
        splice.at[name, "reporter_p_corrected"] if name in splice.index else 1.0))[:12][::-1]
    y = np.arange(len(ranked))
    fig, ax = plt.subplots(figsize=(12.8, max(5.2, 0.48 * len(ranked) + 1.9)))
    for shift, table, color, label in ((-0.18, baseline, TEAL, "Baseline AD vs controls"),
                                       (0.18, splice, ORANGE, "Splicing AD vs controls")):
        values = [float(table.at[name, "reporter_p_corrected"]) if name in table.index else np.nan for name in ranked]
        width = [-np.log10(max(value, 1e-300)) if np.isfinite(value) else np.nan for value in values]
        ax.barh(y + shift, width, height=0.33, color=color, label=label)
    ax.axvline(-np.log10(0.05), color=NAVY, linestyle="--", linewidth=1.2, label="Corrected p = 0.05")
    labels = [name if len(name) <= 72 else name[:69] + "…" for name in ranked]
    ax.set_yticks(y, labels)
    ax.set_xlabel("−log10(reporter-like corrected p-value)")
    ax.set_title("Reporter-like pathway comparison", loc="left", color=NAVY, fontsize=16, pad=14)
    fig.legend(*ax.get_legend_handles_labels(), frameon=False, fontsize=9,
               loc="upper left", bbox_to_anchor=(0.17, 0.92), ncol=3)
    finish(fig, ax, out / "reporter_pathways.svg",
           "Only pathways with corrected p < 0.05 in either computed contrast · archived permutation correction, not BH",
           top=0.82)


def main() -> None:
    style()
    groups, result = load_results()
    targets, target_metrics = target_expression(groups)
    paired, paired_summary = paired_metrics(groups, result)
    associations = fisher_subsystem_table(result)
    summary = {"synthetic": True, "samples": len(groups), "ad_labelled": 100, "control_labelled": 100,
               "splicing_targets": target_metrics, "paired_model_size_and_similarity": paired_summary,
               "conditions": {condition: {
                   "median_active_reactions_ad": float(result[condition]["summary"].loc[
                       result[condition]["summary"].group.eq("ad"), "active_reactions"].median()),
                   "fisher_raw_p_lt_0_05": int(result[condition]["fisher"].p_value.lt(0.05).sum()),
                   "fisher_bh_p_lt_0_05": int(result[condition]["fisher"].p_adjusted_bh.lt(0.05).sum()),
                   "reporter_corrected_p_lt_0_05": int(result[condition]["reporter"].reporter_p_corrected.lt(0.05).sum())}
                   for condition in CONDITIONS}}
    comparison = RESULTS / "comparison"
    comparison.mkdir(parents=True, exist_ok=True)
    targets.to_csv(comparison / "target_expression_categories.tsv", sep="\t", index=False)
    paired.to_csv(comparison / "paired_reaction_metrics.tsv", sep="\t", index=False)
    associations.to_csv(comparison / "fisher_subsystem_associations.tsv", sep="\t", index=False)
    (comparison / "figure_metrics.json").write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    figures = RESULTS / "figures"
    figures.mkdir(parents=True, exist_ok=True)
    plot_activity(groups, result, figures)
    plot_paired_size(paired, paired_summary, figures)
    plot_target_expression(targets, target_metrics, figures)
    plot_reporter(result, figures)
    print(f"Four data-derived SVG figures and audit metrics written to {RESULTS}")


if __name__ == "__main__":
    main()
