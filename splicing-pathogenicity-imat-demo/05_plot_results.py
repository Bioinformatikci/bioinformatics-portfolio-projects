"""Create compact figures from the synthetic demonstration result tables."""
from __future__ import annotations

import matplotlib.pyplot as plt
import pandas as pd

from pipeline_common import CONDITIONS, RESULTS, require


DISPLAY = {
    "expression_only": "Expression only",
    "splicing_informed": "Splicing informed",
    "splicing_pathogenicity_informed": "Splicing + pathogenicity informed",
}


def main() -> None:
    figures = RESULTS / "figures"; figures.mkdir(parents=True, exist_ok=True)
    summaries, active_rows = [], []
    reporters = []
    for condition in CONDITIONS:
        matrix_summary = pd.read_csv(require(RESULTS / condition / "matrix" / "sample_summary.tsv"), sep="\t")
        fisher = pd.read_csv(require(RESULTS / condition / "fisher" / "fisher_all_reactions.tsv"), sep="\t")
        reporter = pd.read_csv(require(RESULTS / condition / "reporter_like" / "reporter_like_all_pathways.tsv"), sep="\t")
        active_rows.extend({"condition": DISPLAY[condition], "group": row.group, "active_reactions": row.active_reactions}
                           for row in matrix_summary.itertuples(index=False))
        summaries.append({"condition": DISPLAY[condition], "fisher_p_lt_0_05": int(fisher.p_value.lt(0.05).sum()),
                          "reporter_like_adjusted_p_lt_0_05": int(reporter.reporter_p_adjusted_bh.lt(0.05).sum())})
        reporters.append(reporter.assign(condition=DISPLAY[condition]))

    active = pd.DataFrame(active_rows)
    fig, axis = plt.subplots(figsize=(9, 4.5))
    positions, labels = [], []
    for i, condition in enumerate(DISPLAY.values()):
        for j, group in enumerate(("control", "ad")):
            positions.append(i * 3 + j + 1); labels.append(f"{condition}\n{group}")
    values = [active.loc[(active.condition.eq(condition)) & (active.group.eq(group)), "active_reactions"]
              for condition in DISPLAY.values() for group in ("control", "ad")]
    axis.boxplot(values, positions=positions)
    axis.set(xticks=positions, xticklabels=labels, ylabel="Active reaction count")
    axis.tick_params(axis="x", rotation=18)
    fig.tight_layout(); fig.savefig(figures / "active_reaction_counts.png", dpi=220); plt.close(fig)

    counts = pd.DataFrame(summaries).set_index("condition")
    axis = counts.plot.bar(figsize=(8, 4.5), color=["#365f91", "#bd7b32"])
    axis.set(ylabel="Count", xlabel="", title="Significant result counts")
    axis.legend(["Fisher p < 0.05", "Reporter-like adjusted p < 0.05"])
    axis.figure.tight_layout(); axis.figure.savefig(figures / "significant_result_counts.png", dpi=220); plt.close(axis.figure)

    reporter_all = pd.concat(reporters, ignore_index=True).sort_values("reporter_p_adjusted_bh").head(12)
    if reporter_all.empty:
        raise RuntimeError("No reporter-like pathways were available to plot.")
    fig, axis = plt.subplots(figsize=(9, 5))
    axis.barh(reporter_all.subsystem, -__import__("numpy").log10(reporter_all.reporter_p_adjusted_bh.clip(lower=1e-300)), color="#587a5b")
    axis.set(xlabel="-log10(BH-adjusted reporter-like p-value)", title="Top reporter-like pathways (synthetic demo)")
    axis.invert_yaxis(); fig.tight_layout(); fig.savefig(figures / "top_reporter_like_pathways.png", dpi=220); plt.close(fig)
    print(f"Figures written to {figures}")


if __name__ == "__main__":
    main()
