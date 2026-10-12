"""Run Fisher and corrected reporter-like analysis for completed layers."""
from __future__ import annotations

import argparse
import html
import subprocess

import numpy as np

import pandas as pd

from pipeline_common import (CONDITIONS, RAW, RESULTS, ROOT, SEED, _bh, condition_dir,
                             require, run_fisher)


def direction_table(matrix: pd.DataFrame, groups: pd.DataFrame) -> pd.Series:
    ad = groups.loc[groups.group.eq("ad"), "sample_id"].astype(str)
    control = groups.loc[groups.group.eq("control"), "sample_id"].astype(str)
    return matrix[ad].mean(axis=1).sub(matrix[control].mean(axis=1)).rename("ad_minus_control_activity")


def compare_fisher(left: pd.DataFrame, right: pd.DataFrame, left_name: str, right_name: str) -> pd.DataFrame:
    merged = left.merge(right, on="reaction_id", suffixes=(f"_{left_name}", f"_{right_name}"), validate="one_to_one")
    left_sig = merged[f"p_value_{left_name}"].lt(0.05)
    right_sig = merged[f"p_value_{right_name}"].lt(0.05)
    merged["significance_membership"] = "neither"
    merged.loc[left_sig & right_sig, "significance_membership"] = "both"
    merged.loc[left_sig & ~right_sig, "significance_membership"] = left_name + "_only"
    merged.loc[~left_sig & right_sig, "significance_membership"] = right_name + "_only"
    return merged


def compare_reporter(left: pd.DataFrame, right: pd.DataFrame, left_name: str, right_name: str) -> pd.DataFrame:
    merged = left.merge(right, on="subsystem", how="outer",
                        suffixes=(f"_{left_name}", f"_{right_name}"), validate="one_to_one")
    left_sig = merged[f"reporter_p_corrected_{left_name}"].lt(0.05)
    right_sig = merged[f"reporter_p_corrected_{right_name}"].lt(0.05)
    merged["significance_membership"] = "neither"
    merged.loc[left_sig & right_sig, "significance_membership"] = "both"
    merged.loc[left_sig & ~right_sig, "significance_membership"] = left_name + "_only"
    merged.loc[~left_sig & right_sig, "significance_membership"] = right_name + "_only"
    return merged


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--through", choices=CONDITIONS, default=CONDITIONS[-1],
                        help="Stop after this completed layer; the default analyses all three")
    args = parser.parse_args()
    selected = CONDITIONS[:CONDITIONS.index(args.through) + 1]
    subsystem = pd.read_csv(require(RAW / "reaction_subsystem.tsv"), sep="\t")
    if subsystem.columns.tolist() != ["reaction_id", "SUBSYSTEM"] or len(subsystem) != 13071:
        raise ValueError("Human-GEM subsystem annotation differs from the archived workflow")
    condition_data, summary_rows = {}, []
    for condition in selected:
        output = condition_dir(condition)
        matrix = pd.read_csv(require(output / "matrix" / "reaction_activity_matrix.tsv"), sep="\t", index_col=0)
        groups = pd.read_csv(require(RESULTS / "expression_only" / "processed" / "sample_groups.tsv"), sep="\t")
        fisher_dir, reporter_dir = output / "fisher", output / "reporter_like"
        fisher_dir.mkdir(exist_ok=True); reporter_dir.mkdir(exist_ok=True)
        fisher = run_fisher(matrix, groups)
        fisher["ad_minus_control_activity"] = direction_table(matrix, groups).reindex(fisher.reaction_id).to_numpy()
        fisher.to_csv(fisher_dir / "fisher_all_reactions.tsv", sep="\t", index=False)
        fisher.loc[fisher.p_value.lt(0.05)].to_csv(fisher_dir / "fisher_p_lt_0_05.tsv", sep="\t", index=False)
        fisher.loc[fisher.p_adjusted_bh.lt(0.05)].to_csv(fisher_dir / "fisher_adjusted_p_lt_0_05.tsv", sep="\t", index=False)
        source_columns = ["reaction_id", "ad_present", "ad_absent", "control_present",
                          "control_absent", "odds_ratio", "p_value"]
        fisher[source_columns].to_csv(fisher_dir / "fisher_for_reporter.tsv", sep="\t", index=False)
        subprocess.run(["Rscript", str(ROOT / "reporter_like.R"),
                        str(fisher_dir / "fisher_for_reporter.tsv"),
                        str(RAW / "reaction_subsystem.tsv"), str(reporter_dir), str(SEED)], check=True)
        original_reporter = pd.read_csv(reporter_dir / "reporter_like.tsv", sep="\t")
        reporter = original_reporter.rename(columns={
            "SUBSYSTEM": "subsystem", "freq": "reaction_count", "zsum": "z_sum",
            "p_path": "reporter_score", "uncorrected_pvalue_path": "reporter_p_uncorrected",
            "mean_correction": "permutation_mean", "sd_correction": "permutation_sd",
            "p_path_corrected": "reporter_score_corrected",
            "pvalue_path_corrected": "reporter_p_corrected",
        })
        for column in reporter.columns.difference(["subsystem"]):
            reporter[column] = pd.to_numeric(reporter[column], errors="coerce")
        reporter["reporter_p_adjusted_bh"] = np.nan
        finite = reporter.reporter_p_corrected.notna()
        if finite.any():
            reporter.loc[finite, "reporter_p_adjusted_bh"] = _bh(
                reporter.loc[finite, "reporter_p_corrected"].to_numpy(dtype=float))
        reporter.to_csv(reporter_dir / "reporter_like_all_pathways.tsv", sep="\t", index=False)
        reporter.loc[reporter.reporter_p_adjusted_bh.lt(0.05)].to_csv(reporter_dir / "reporter_like_adjusted_p_lt_0_05.tsv", sep="\t", index=False)
        condition_data[condition] = {"fisher": fisher, "reporter": reporter, "matrix": matrix}
        summary_rows.append({"condition": condition, "reactions": len(matrix),
                             "fisher_raw_p_lt_0_05": int(fisher.p_value.lt(0.05).sum()),
                             "fisher_bh_p_lt_0_05": int(fisher.p_adjusted_bh.lt(0.05).sum()),
                             "reporter_like_corrected_p_lt_0_05": int(reporter.reporter_p_corrected.lt(0.05).sum()),
                             "reporter_like_bh_p_lt_0_05": int(reporter.reporter_p_adjusted_bh.lt(0.05).sum())})

    comparison = RESULTS / "comparison"; comparison.mkdir(parents=True, exist_ok=True)
    pairs = [("expression_only", "splicing_informed"), ("splicing_informed", "splicing_pathogenicity_informed")]
    for left_name, right_name in pairs:
        if left_name not in condition_data or right_name not in condition_data:
            continue
        compare_fisher(condition_data[left_name]["fisher"], condition_data[right_name]["fisher"], left_name, right_name).to_csv(
            comparison / f"fisher_{left_name}_vs_{right_name}.tsv", sep="\t", index=False)
        compare_reporter(condition_data[left_name]["reporter"], condition_data[right_name]["reporter"], left_name, right_name).to_csv(
            comparison / f"reporter_like_{left_name}_vs_{right_name}.tsv", sep="\t", index=False)
    summary = pd.DataFrame(summary_rows)
    complete = len(selected) == len(CONDITIONS)
    summary.to_csv(comparison / ("analysis_summary.tsv" if complete else "analysis_summary_partial.tsv"),
                   sep="\t", index=False)
    if complete:
        write_html_report(summary, comparison)
    print(f"Fisher and reporter-like results for {len(selected)} completed layers written to {RESULTS}")


def write_html_report(summary: pd.DataFrame, comparison: pd.PathLike) -> None:
    rows = "\n".join("<tr>" + "".join(f"<td>{html.escape(str(value))}</td>" for value in row) + "</tr>"
                     for row in summary.itertuples(index=False, name=None))
    headers = "".join(f"<th>{html.escape(column)}</th>" for column in summary.columns)
    report = f"""<!doctype html><html lang=\"en\"><head><meta charset=\"utf-8\"><title>One-cohort iMAT demo report</title>
<style>body{{font:16px system-ui,sans-serif;line-height:1.5;max-width:900px;margin:40px auto;padding:0 20px;color:#16212d}}table{{border-collapse:collapse;width:100%}}th,td{{padding:8px;border:1px solid #ccd5df;text-align:left}}th{{background:#eaf1f8}}code{{background:#f0f3f5;padding:2px 4px}}</style></head><body>
<h1>One-cohort splicing and pathogenicity iMAT demonstration</h1>
<p>This report describes a synthetic-data demonstration; it contains no study cohort results. Reaction scores use COBAMP GPR aggregation and binary reaction matrices use TROPPO iMAT. Context-specific GEMs are not exported.</p>
<ol><li>Expression-only reaction activity was inferred.</li><li>Only AD splicing-target expression values were set to zero and quality-controlled.</li><li>Only AD genes with <code>pathogenicity_rank &gt; 0.95</code> were additionally set to zero and quality-controlled.</li><li>All reaction-level AD/control Fisher tests and permutation reporter-like pathway tests were saved.</li></ol>
<p>Reporter-like analysis follows the archived R implementation: raw Fisher p &lt; 0.05 reactions, 1,000 within-table permutations, seed 20261010, and a permutation-normalized pathway score. BH adjustments are additional columns for transparency; they are not part of the archived reporter calculation. Hypergeometric analysis was not performed.</p>
<h2>Result counts</h2><table><thead><tr>{headers}</tr></thead><tbody>{rows}</tbody></table>
<p>Detailed result tables are in the sibling condition folders; condition comparisons are in this folder.</p></body></html>"""
    (comparison / "analysis_report.html").write_text(report, encoding="utf-8")


if __name__ == "__main__":
    main()
