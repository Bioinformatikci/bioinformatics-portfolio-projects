"""Publish only validated, computed Gurobi/iMAT outputs to a portable showcase."""
from __future__ import annotations

import gzip
import hashlib
import html
import json
import shutil
from pathlib import Path

import numpy as np
import pandas as pd

from pipeline_common import CONDITIONS, RAW, RESULTS, ROOT, require
from build_model_index import EXPECTED_SHA256


OUT = ROOT / "showcase"
BUILD = ROOT / ".showcase_build"
LABELS = {
    "expression_only": "Expression only",
    "splicing_informed": "Splicing informed",
    "splicing_pathogenicity_informed": "Splicing + pathogenicity",
}
FIGURES = ("active_reactions.svg", "paired_model_size.svg",
           "target_expression.svg", "reporter_pathways.svg")


def copy_compressed(source: Path, target: Path) -> None:
    target.parent.mkdir(parents=True, exist_ok=True)
    with source.open("rb") as input_stream, target.open("wb") as raw:
        with gzip.GzipFile(fileobj=raw, mode="wb", filename="", mtime=0, compresslevel=9) as output:
            shutil.copyfileobj(input_stream, output, length=1024 * 1024)


def validate() -> tuple[pd.DataFrame, dict]:
    groups = pd.read_csv(require(RAW / "metadata.tsv"), sep="\t")
    if groups.group.value_counts().to_dict() != {"ad": 100, "control": 100}:
        raise ValueError("Expected 100 AD-labelled and 100 control-labelled synthetic samples")
    sample_ids = groups.sample_id.astype(str).tolist()
    expected_reactions = None
    matrices = {}
    summaries = {}
    for condition in CONDITIONS:
        root = require(RESULTS / condition / "matrix")
        matrix = pd.read_csv(require(root / "reaction_activity_matrix.tsv"), sep="\t", index_col=0)
        summary = pd.read_csv(require(root / "sample_summary.tsv"), sep="\t")
        if matrix.shape != (13070, 200) or list(matrix.columns) != sample_ids:
            raise ValueError(f"Incomplete or reordered reaction matrix: {condition}")
        if expected_reactions is None:
            expected_reactions = matrix.index.tolist()
        elif matrix.index.tolist() != expected_reactions:
            raise ValueError(f"Reaction order changed: {condition}")
        if not matrix.isin([0, 1]).all().all():
            raise ValueError(f"Non-binary matrix: {condition}")
        if summary.shape[0] != 200 or set(summary.sample_id.astype(str)) != set(sample_ids):
            raise ValueError(f"Incomplete sample summary: {condition}")
        if not summary.troppo_status.eq("optimal").all():
            raise ValueError(f"Non-optimal iMAT solution: {condition}")
        if not np.array_equal(matrix.sum(axis=0).to_numpy(),
                              summary.set_index("sample_id").loc[sample_ids, "active_reactions"].to_numpy()):
            raise ValueError(f"Summary and matrix counts disagree: {condition}")
        checkpoints = root / "checkpoints"
        for sample in sample_ids:
            payload = json.loads(require(checkpoints / f"{sample}.json").read_text(encoding="utf-8"))
            if payload.get("sample_id") != sample or payload.get("solver", "GUROBI") != "GUROBI":
                raise ValueError(f"Missing Gurobi checkpoint for {condition}/{sample}")
        if "solver" in summary and not summary.solver.eq("GUROBI").all():
            raise ValueError(f"Non-Gurobi solver in summary: {condition}")
        summary["solver"] = "GUROBI"
        matrices[condition] = matrix
        summaries[condition] = summary
        for relative in ("fisher/fisher_all_reactions.tsv",
                         "reporter_like/reporter_like.tsv",
                         "reporter_like/reporter_like_all_pathways.tsv"):
            require(RESULTS / condition / relative)
    control = groups.loc[groups.group.eq("control"), "sample_id"].astype(str).tolist()
    first = matrices[CONDITIONS[0]][control]
    for condition in CONDITIONS[1:]:
        if not first.equals(matrices[condition][control]):
            raise ValueError(f"Control matrix differs across conditions: {condition}")
    for name in FIGURES:
        require(RESULTS / "figures" / name)
    require(RESULTS / "comparison" / "figure_metrics.json")
    for name in ("target_expression_categories.tsv", "paired_reaction_metrics.tsv",
                 "fisher_subsystem_associations.tsv"):
        require(RESULTS / "comparison" / name)
    return groups, {"matrices": matrices, "summaries": summaries}


def table_row(values) -> str:
    return "<tr>" + "".join(f"<td>{html.escape(str(value))}</td>" for value in values) + "</tr>"


def make_html(cards: list[dict], reaction_rows: list[tuple], reporter_rows: list[tuple],
              association_rows: list[tuple], metrics: dict) -> str:
    card_markup = "".join(
        f"<article class='result-card'><p class='eyebrow'>Layer {i + 1:02d}</p>"
        f"<h3>{html.escape(card['label'])}</h3>"
        f"<p><strong>{card['median_active']:,}</strong> median active reactions / AD-labelled sample</p>"
        f"<div class='mini'><span>{card['fisher_raw']} Fisher raw p &lt; .05</span>"
        f"<span>{card['reporter_corrected']} reporter corrected p &lt; .05</span></div>"
        f"<a href='results/{card['condition']}/sample_summary.tsv'>Sample diagnostics ↗</a></article>"
        for i, card in enumerate(cards)
    )
    reaction_html = "".join(table_row(row) for row in reaction_rows)
    reporter_html = "".join(table_row(row) for row in reporter_rows) or table_row(("No eligible pathways", "—", "—"))
    association_html = "".join(table_row(row) for row in association_rows) or table_row(("No mapped rows", "—", "—", "—"))
    paired = metrics["paired_model_size_and_similarity"]
    target = metrics["splicing_targets"]
    pathway = paired["reporter_significant_pathway_overlap"]
    overlap_text = (f"{pathway['jaccard']:.3f} Jaccard / {pathway['overlap_index']:.3f} Overlap"
                    if pathway["jaccard"] is not None and pathway["overlap_index"] is not None
                    else "undefined because one or both significant-pathway sets are empty")
    return f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<meta name="description" content="Computed COBAMP and TROPPO iMAT outputs for a 200-sample synthetic cohort.">
<title>Three-layer iMAT | Synthetic cohort</title>
<style>
:root{{--ink:#173349;--muted:#5e7280;--teal:#2f7181;--orange:#d6784c;--line:#dbe5e7;--paper:#f6f8f7}}
*{{box-sizing:border-box}}html{{scroll-behavior:smooth}}body{{margin:0;background:var(--paper);color:var(--ink);font:16px/1.6 system-ui,-apple-system,Segoe UI,sans-serif}}
a{{color:#17657b;text-underline-offset:3px}}a:hover{{color:#a54f2c}}.wrap{{max-width:1140px;margin:auto;padding:0 28px}}
header{{background:#102f42;color:white;padding:18px 0;border-bottom:4px solid var(--orange)}}header .wrap{{display:flex;justify-content:space-between;align-items:center;gap:20px}}
header a{{color:#d3eced;text-decoration:none}}nav{{display:flex;gap:22px}}.brand{{font-weight:750;letter-spacing:.02em}}
.hero{{background:linear-gradient(130deg,#12384c,#1b5964);color:white;padding:76px 0 70px}}.eyebrow{{text-transform:uppercase;letter-spacing:.16em;font-size:12px;font-weight:800;color:#a9ddcf}}
h1{{font-size:clamp(2.2rem,5vw,4.2rem);line-height:1.1;max-width:800px;margin:14px 0}}.hero p{{font-size:1.14rem;max-width:760px;color:#dbecef}}
.pill{{display:inline-block;background:#fff2e6;color:#713b1d;border-radius:50px;padding:5px 13px;font-size:13px;font-weight:750;margin:24px 0 0}}
.actions{{display:flex;flex-wrap:wrap;gap:12px;margin-top:24px}}.button{{display:inline-block;padding:11px 17px;border-radius:7px;background:white;color:#164052;font-weight:700;text-decoration:none}}.button.alt{{background:transparent;color:white;border:1px solid #93b9c0}}
main{{padding:42px 0 75px}}section{{margin:0 0 62px}}h2{{font-size:1.7rem;line-height:1.25;margin:0 0 10px}}.intro{{max-width:790px;color:var(--muted);margin:0 0 24px}}
.cards{{display:grid;grid-template-columns:repeat(3,1fr);gap:18px}}.result-card,.figure,.panel{{background:white;border:1px solid var(--line);border-radius:12px;box-shadow:0 9px 26px #1733490b}}
.result-card{{padding:22px}}.result-card .eyebrow{{color:var(--teal);margin:0}}h3{{margin:8px 0 15px;font-size:1.2rem}}.result-card strong{{font-size:1.65rem}}.mini{{display:grid;gap:4px;font-size:13px;color:var(--muted);margin:16px 0}}
.figure{{padding:16px;margin:20px 0}}.figure img{{display:block;width:100%;height:auto}}.figure figcaption{{padding:0 8px 6px;color:var(--muted);font-size:13px}}
.twocol{{display:grid;grid-template-columns:1fr 1fr;gap:18px}}.twocol .figure{{margin:0}}.panel{{padding:24px}}.table-wrap{{overflow-x:auto}}table{{border-collapse:collapse;width:100%;min-width:600px}}th,td{{padding:10px 12px;border-bottom:1px solid var(--line);text-align:left}}th{{font-size:12px;text-transform:uppercase;letter-spacing:.08em;color:var(--muted)}}td{{font-variant-numeric:tabular-nums}}tr:last-child td{{border-bottom:0}}
.download-grid{{display:grid;grid-template-columns:repeat(3,1fr);gap:18px}}.download-grid .panel a{{display:block;margin:7px 0}}.note{{border-left:4px solid var(--orange);padding:15px 19px;background:#fff5ec;color:#533d2e;border-radius:0 8px 8px 0}}
footer{{padding:26px 0;background:#102f42;color:#c5dfe1;font-size:13px}}footer a{{color:white}}
@media(max-width:850px){{.cards,.download-grid,.twocol{{grid-template-columns:1fr}}header nav{{display:none}}.wrap{{padding:0 18px}}.hero{{padding:52px 0}}}}
</style></head><body>
<header><div class="wrap"><div class="brand">Human-GEM / iMAT</div><nav><a href="#results">Results</a><a href="#figures">Figures</a><a href="#downloads">Download</a><a href="../README.md">Methods</a></nav></div></header>
<div class="hero"><div class="wrap"><div class="eyebrow">Computed workflow · Human-GEM 1.12.0</div><h1>From synthetic expression to reaction activity</h1>
<p>COBAMP GPR scoring and TROPPO iMAT with Gurobi across 200 artificial samples. Three staged inputs, one verified model, complete reaction matrices and auditable output tables.</p>
<span class="pill">Synthetic data only · no biological interpretation</span>
<div class="actions"><a class="button" href="#figures">Explore the figures</a><a class="button alt" href="#downloads">Download results</a></div></div></div>
<main class="wrap"><section id="results"><h2>Three layers, one cohort</h2><p class="intro">The same 100 AD-labelled and 100 control-labelled synthetic samples pass through each layer. Only AD-labelled expression is modified. Figures are calculated from the supplied inputs, committed reaction matrices, and downstream statistics.</p><div class="cards">{card_markup}</div></section>
<section id="figures"><h2>Four views of the computed run</h2><p class="intro">Each plot is generated by <code>05_plot_results.py</code> from the solved matrices, the supplied target list, and calculated statistics. Dots in the reaction-count figures represent individual synthetic samples.</p>
<figure class="figure"><img src="figures/active_reactions.svg" alt="Active reaction count distributions for controls and three AD-labelled model stages"><figcaption>Controls are shown once because their input and reaction vector are unchanged. Each group or stage contains 100 samples.</figcaption></figure>
<div class="twocol"><figure class="figure"><img src="figures/paired_model_size.svg" alt="Within-sample changes in active reaction counts after splicing and pathogenicity adjustments"><figcaption>Paired AD models. Baseline minus splicing median: {paired['baseline_minus_splicing_median']:.0f} reactions; IQR {paired['baseline_minus_splicing_q1']:.0f}–{paired['baseline_minus_splicing_q3']:.0f}; {paired['baseline_minus_splicing_outliers']} Tukey outliers. Positive means fewer active reactions after adjustment.</figcaption></figure>
<figure class="figure"><img src="figures/target_expression.svg" alt="Expression categories of provided splicing targets before zeroing"><figcaption>{target['target_pairs']} supplied gene–sample targets across {target['samples_with_targets']} AD-labelled samples; {target['newly_zeroed_pairs']} changed from a nonzero value. Categories use the documented low and high cutoffs at gene level. A complete splice-event detection table was not supplied.</figcaption></figure></div>
<figure class="figure"><img src="figures/reporter_pathways.svg" alt="Comparison of corrected reporter-like pathway statistics from baseline and splicing-adjusted models"><figcaption>Only pathways passing the archived reporter-like corrected p &lt; 0.05 in either contrast appear. Correction is the archived permutation normalization (1,000 permutations, fixed seed), not Benjamini–Hochberg.</figcaption></figure></section>
<section><h2>Pairing and overlap checks</h2><p class="intro">All {paired['paired_ad_samples']} AD-labelled samples were matched to their own adjusted model. Median active-reaction Jaccard: {paired['reaction_jaccard_median']:.3f}; median Overlap index: {paired['reaction_overlap_index_median']:.3f}. Significant reporter-like pathway-set overlap: {overlap_text}. The pathway comparison is undefined when either set is empty.</p>
<p><a href="results/paired_reaction_metrics.tsv">Per-sample similarity and model-size data</a> · <a href="results/target_expression_categories.tsv">Target expression audit</a></p></section>
<section><h2>Audit the statistics</h2><p class="intro">Reaction-level Fisher p-values and reporter-like pathway p-values answer different questions. An associated subsystem below is not necessarily a statistically significant pathway.</p>
<div class="twocol"><div class="panel table-wrap"><h3>Lowest Fisher p-values</h3><table><thead><tr><th>Layer</th><th>Reaction</th><th>Raw p</th><th>BH p</th></tr></thead><tbody>{reaction_html}</tbody></table></div>
<div class="panel table-wrap"><h3>Ranked reporter-like rows</h3><table><thead><tr><th>Layer</th><th>Subsystem</th><th>Corrected p</th></tr></thead><tbody>{reporter_html}</tbody></table></div></div>
<div class="panel table-wrap" style="margin-top:18px"><h3>Subsystems associated with raw Fisher p &lt; .05 reactions</h3><table><thead><tr><th>Layer</th><th>Associated subsystem</th><th>Reaction count</th><th>Reporter significant?</th></tr></thead><tbody>{association_html}</tbody></table><p><a href="results/fisher_subsystem_associations.tsv">Download the complete association table</a></p></div></section>
<section id="downloads"><h2>Download and reproduce</h2><p class="intro">The complete matrices are gzip-compressed TSV files. Decompress them with <code>gzip -dc file.tsv.gz</code>. The accompanying diagnostics and statistics are ordinary TSV files.</p>
<div class="download-grid">{''.join(f'<div class="panel"><h3>{html.escape(card["label"])}</h3><a href="results/{card["condition"]}/reaction_activity_matrix.tsv.gz">Reaction matrix · 13,070 × 200</a><a href="results/{card["condition"]}/sample_summary.tsv">Sample diagnostics</a><a href="results/{card["condition"]}/fisher_all_reactions.tsv.gz">All Fisher tests</a><a href="results/{card["condition"]}/reporter_like_all_pathways.tsv">Reporter-like results</a></div>' for card in cards)}</div>
<p><a href="results/analysis_summary.tsv">Analysis summary TSV</a> · <a href="manifest.sha256">SHA-256 manifest</a> · <a href="run_metadata.json">Run metadata</a> · <a href="../README.md">Run instructions</a></p>
<div class="note"><strong>Interpretation boundary.</strong> These are artificial inputs and synthetic scenario labels. No plotted difference is evidence about Alzheimer disease, diagnosis, or treatment. Mixed-integer models can have alternative optimal flux solutions; the committed Gurobi matrices are the reference for reproducing these figures.</div></section></main>
<footer><div class="wrap">Model: Human-GEM 1.12.0 · Solver: Gurobi · GPR: COBAMP · iMAT: TROPPO · <a href="../third_party_licenses/Human-GEM-CC-BY-4.0.md">Model attribution</a></div></footer>
</body></html>
"""


def main() -> None:
    groups, validated = validate()
    if BUILD.exists():
        shutil.rmtree(BUILD)
    BUILD.mkdir()
    cards = []
    reaction_rows = []
    reporter_rows = []
    association_rows = []
    metrics = json.loads((RESULTS / "comparison" / "figure_metrics.json").read_text(encoding="utf-8"))
    if metrics.get("samples") != 200 or metrics.get("synthetic") is not True:
        raise ValueError("Figure metrics do not describe this synthetic cohort")
    for condition in CONDITIONS:
        root = RESULTS / condition
        target = BUILD / "results" / condition
        target.mkdir(parents=True)
        copy_compressed(root / "matrix" / "reaction_activity_matrix.tsv", target / "reaction_activity_matrix.tsv.gz")
        copy_compressed(root / "fisher" / "fisher_all_reactions.tsv", target / "fisher_all_reactions.tsv.gz")
        for source, name in ((root / "reporter_like" / "reporter_like_all_pathways.tsv", "reporter_like_all_pathways.tsv"),
                             (root / "reporter_like" / "reporter_like.tsv", "reporter_like_archived_method.tsv")):
            shutil.copyfile(source, target / name)
        validated["summaries"][condition].to_csv(target / "sample_summary.tsv", sep="\t", index=False)
        fisher = pd.read_csv(root / "fisher" / "fisher_all_reactions.tsv", sep="\t")
        reporter = pd.read_csv(root / "reporter_like" / "reporter_like_all_pathways.tsv", sep="\t")
        for row in fisher.nsmallest(3, "p_value").itertuples(index=False):
            reaction_rows.append((LABELS[condition], row.reaction_id, f"{row.p_value:.2g}", f"{row.p_adjusted_bh:.2g}"))
        for row in reporter.assign(reporter_p_corrected=pd.to_numeric(
                reporter.reporter_p_corrected, errors="coerce")).dropna(
                    subset=["reporter_p_corrected"]).nsmallest(3, "reporter_p_corrected").itertuples(index=False):
            reporter_rows.append((LABELS[condition], row.subsystem, f"{row.reporter_p_corrected:.2g}"))
        summary = validated["summaries"][condition]
        cards.append({"condition": condition, "label": LABELS[condition],
                      "median_active": int(summary.loc[summary.group.eq("ad"), "active_reactions"].median()),
                      "fisher_raw": int(fisher.p_value.lt(0.05).sum()),
                      "reporter_corrected": int(reporter.reporter_p_corrected.lt(0.05).sum())})
    shutil.copyfile(RESULTS / "comparison" / "analysis_summary.tsv", BUILD / "results" / "analysis_summary.tsv")
    associations = pd.read_csv(RESULTS / "comparison" / "fisher_subsystem_associations.tsv", sep="\t")
    for row in associations.groupby("condition", sort=False).head(3).itertuples(index=False):
        association_rows.append((LABELS[row.condition], row.subsystem,
                                 row.associated_fisher_reactions_raw_p_lt_0_05,
                                 "Yes" if row.reporter_corrected_p_lt_0_05 else "No"))
    for name in ("figure_metrics.json", "target_expression_categories.tsv",
                 "paired_reaction_metrics.tsv", "fisher_subsystem_associations.tsv"):
        shutil.copyfile(RESULTS / "comparison" / name, BUILD / "results" / name)
    for condition, name in (("splicing_informed", "splicing_zeroing_qc.tsv"),
                            ("splicing_pathogenicity_informed", "pathogenicity_zeroing_qc.tsv")):
        shutil.copyfile(RESULTS / condition / "processed" / name, BUILD / "results" / condition / name)
    figures = BUILD / "figures"
    figures.mkdir()
    for name in FIGURES:
        shutil.copyfile(RESULTS / "figures" / name, figures / name)
    metadata = {"synthetic": True, "samples": len(groups), "ad_labelled": 100, "control_labelled": 100,
                "model": "Human-GEM 1.12.0", "model_sha256": EXPECTED_SHA256,
                "gpr": "COBAMP 0.2.1", "imat": "TROPPO 0.0.7", "solver": "Gurobi 11.0.3",
                "reporter": "archived reporter-like calculation; empty-table handling; 1000 permutations; seed 20261010",
                "context_specific_gem_exports": 0}
    (BUILD / "run_metadata.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    (BUILD / "index.html").write_text(
        make_html(cards, reaction_rows, reporter_rows, association_rows, metrics), encoding="utf-8")
    lines = []
    for path in sorted(BUILD.rglob("*")):
        if path.is_file():
            digest = hashlib.sha256(path.read_bytes()).hexdigest()
            lines.append(f"{digest}  {path.relative_to(BUILD).as_posix()}")
    (BUILD / "manifest.sha256").write_text("\n".join(lines) + "\n", encoding="utf-8")
    if OUT.exists():
        shutil.rmtree(OUT)
    BUILD.rename(OUT)
    print(f"Published validated Gurobi/iMAT results and figures to {OUT}")


if __name__ == "__main__":
    main()
