"""Run Fisher and corrected reporter-like analysis for all three layers."""
from __future__ import annotations

import html

import pandas as pd

from pipeline_common import (CONDITIONS, RAW, RESULTS, condition_dir, reporter_like, require,
                             run_fisher)


def direction_table(matrix: pd.DataFrame, groups: pd.DataFrame) -> pd.Series:
    ad = groups.loc[groups.group.eq("ad"), "sample_id"].astype(str)
    control = groups.loc[groups.group.eq("control"), "sample_id"].astype(str)
    return matrix[ad].mean(axis=1).sub(matrix[control].mean(axis=1)).rename("ad_minus_control_activity")


def compare_fisher(left: pd.DataFrame, right: pd.DataFrame, left_name: str, right_name: str) -> pd.DataFrame:
    merged = left.merge(right, on="reaction_id", suffixes=(f"_{left_name}", f"_{right_name}"), validate="one_to_one")
    left_sig, right_sig = merged[f"p_value_{left_name}"].lt(0.05), merged[f"p_value_{right_name}"].lt(0.05)
    merged["significance_membership"] = "neither"
    merged.loc[left_sig & right_sig, "significance_membership"] = "both"
    merged.loc[left_sig & ~right_sig, "significance_membership"] = left_name + "_only"
    merged.loc[~left_sig & right_sig, "significance_membership"] = right_name + "_only"
    return merged


def compare_reporter(left: pd.DataFrame, right: pd.DataFrame, left_name: str, right_name: str) -> pd.DataFrame:
    merged = left.merge(right, on="subsystem", suffixes=(f"_{left_name}", f"_{right_name}"), validate="one_to_one")
    left_sig = merged[f"reporter_p_adjusted_bh_{left_name}"].lt(0.05)
    right_sig = merged[f"reporter_p_adjusted_bh_{right_name}"].lt(0.05)
    merged["significance_membership"] = "neither"
    merged.loc[left_sig & right_sig, "significance_membership"] = "both"
    merged.loc[left_sig & ~right_sig, "significance_membership"] = left_name + "_only"
    merged.loc[~left_sig & right_sig, "significance_membership"] = right_name + "_only"
    return merged


def main() -> None:
    subsystem = pd.read_csv(require(RAW / "reaction_subsystem.tsv"), sep="\t")
    condition_data, summary_rows = {}, []
    for condition in CONDITIONS:
        output = condition_dir(condition)
        matrix = pd.read_csv(require(output / "matrix" / "reaction_activity_matrix.tsv"), sep="\t", index_col=0)
        groups = pd.read_csv(require(RESULTS / "expression_only" / "processed" / "sample_groups.tsv"), sep="\t")
        fisher_dir, reporter_dir = output / "fisher", output / "reporter_like"
        fisher_dir.mkdir(exist_ok=True); reporter_dir.mkdir(exist_ok=True)
        fisher = run_fisher(matrix, groups)
        fisher["ad_minus_control_activity"] = direction_table(matrix, groups).reindex(fisher.reaction_id).to_numpy()
        fisher.to_csv(fisher_dir / "fisher_all_reactions.tsv", sep="\t", index=False)
        fisher.loc[fisher.p_value.lt(0.05)].to_csv(fisher_dir / "fisher_p_lt_0_05.tsv", sep="\t", index=False)
        reporter = reporter_like(fisher, subsystem)
        reporter.to_csv(reporter_dir / "reporter_like_all_pathways.tsv", sep="\t", index=False)
        reporter.loc[reporter.reporter_p_adjusted_bh.lt(0.05)].to_csv(reporter_dir / "reporter_like_adjusted_p_lt_0_05.tsv", sep="\t", index=False)
        condition_data[condition] = {"fisher": fisher, "reporter": reporter, "matrix": matrix}
        summary_rows.append({"condition": condition, "reactions": len(matrix),
                             "fisher_p_lt_0_05": int(fisher.p_value.lt(0.05).sum()),
                             "reporter_like_adjusted_p_lt_0_05": int(reporter.reporter_p_adjusted_bh.lt(0.05).sum())})

    comparison = RESULTS / "comparison"; comparison.mkdir(parents=True, exist_ok=True)
    pairs = [("expression_only", "splicing_informed"), ("splicing_informed", "splicing_pathogenicity_informed")]
    for left_name, right_name in pairs:
        compare_fisher(condition_data[left_name]["fisher"], condition_data[right_name]["fisher"], left_name, right_name).to_csv(
            comparison / f"fisher_{left_name}_vs_{right_name}.tsv", sep="\t", index=False)
        compare_reporter(condition_data[left_name]["reporter"], condition_data[right_name]["reporter"], left_name, right_name).to_csv(
            comparison / f"reporter_like_{left_name}_vs_{right_name}.tsv", sep="\t", index=False)
    summary = pd.DataFrame(summary_rows)
    summary.to_csv(comparison / "analysis_summary.tsv", sep="\t", index=False)
    write_html_report(summary, comparison)
    print(f"Fisher and reporter-like results written to {RESULTS}")


def write_html_report(summary: pd.DataFrame, comparison: pd.PathLike) -> None:
    rows = "\n".join("<tr>" + "".join(f"<td>{html.escape(str(value))}</td>" for value in row) + "</tr>"
                     for row in summary.itertuples(index=False, name=None))
    headers = "".join(f"<th>{html.escape(column)}</th>" for column in summary.columns)
    report = f"""<!doctype html><html lang=\"en\"><head><meta charset=\"utf-8\"><title>One-cohort iMAT demo report</title>
<style>body{{font:16px system-ui,sans-serif;line-height:1.5;max-width:900px;margin:40px auto;padding:0 20px;color:#16212d}}table{{border-collapse:collapse;width:100%}}th,td{{padding:8px;border:1px solid #ccd5df;text-align:left}}th{{background:#eaf1f8}}code{{background:#f0f3f5;padding:2px 4px}}</style></head><body>
<h1>One-cohort splicing and pathogenicity iMAT demonstration</h1>
<p>This report describes a synthetic-data demonstration; it contains no study cohort results. Reaction scores use COBAMP GPR aggregation and binary reaction matrices use TROPPO iMAT. Context-specific GEMs are not exported.</p>
<ol><li>Expression-only reaction activity was inferred.</li><li>Only AD splicing-target expression values were set to zero and quality-controlled.</li><li>Only AD genes with <code>pathogenicity_rank &gt; 0.95</code> were additionally set to zero and quality-controlled.</li><li>All reaction-level AD/control Fisher tests and permutation reporter-like pathway tests were saved.</li></ol>
<p>Reporter-like analysis uses 1,000 permutations, seed 20261010 and BH-corrected p-values. Hypergeometric analysis was not performed.</p>
<h2>Result counts</h2><table><thead><tr>{headers}</tr></thead><tbody>{rows}</tbody></table>
<p>Detailed result tables are in the sibling condition folders; condition comparisons are in this folder.</p></body></html>"""
    (comparison / "analysis_report.html").write_text(report, encoding="utf-8")


if __name__ == "__main__":
    main()
