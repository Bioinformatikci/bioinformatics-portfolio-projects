"""Build a deterministic, dependency-free preview of the input/output contract.

This is an illustrative thresholded reaction-state example. It does not run
Human-GEM, COBAMP, TROPPO iMAT, Fisher tests, or reporter-like statistics.
"""
from __future__ import annotations

import argparse
import csv
import html
from pathlib import Path


ROOT = Path(__file__).resolve().parent
SAMPLES = tuple(f"DEMO_AD_{i:02d}" for i in range(1, 5)) + tuple(
    f"DEMO_CONTROL_{i:02d}" for i in range(1, 5)
)
AD = SAMPLES[:4]
CONTROL = SAMPLES[4:]
GENES = tuple(f"SYN_G{i:03d}" for i in range(1, 25))
REACTIONS = tuple(f"SYN_R{i:03d}" for i in range(1, 25))
PATHWAYS = ("Glycolysis", "Lipid metabolism", "Mitochondria", "Amino acids", "Nucleotides", "Transport")
CONDITIONS = ("expression_only", "splicing_informed", "splicing_pathogenicity_informed")
LABELS = ("Expression only", "Splicing informed", "Splicing + pathogenicity")
THRESHOLD = 6.0


def write_tsv(path: Path, headers: tuple[str, ...], rows: list[tuple[object, ...]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(headers)
        writer.writerows(rows)


def svg(path: Path, width: int, height: int, title: str, desc: str, body: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    content = (
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}" role="img" '
        f'aria-labelledby="title desc"><title id="title">{html.escape(title)}</title>'
        f'<desc id="desc">{html.escape(desc)}</desc>'
        f'<rect width="{width}" height="{height}" fill="#ffffff"/>'
        f'<style>text{{font-family:system-ui,Arial,sans-serif;fill:#18332f}}</style>{body}</svg>\n'
    )
    path.write_text(content, encoding="utf-8")


def generate(output: Path) -> None:
    output.mkdir(parents=True, exist_ok=True)
    raw = output / "data" / "raw"
    write_tsv(raw / "metadata.tsv", ("sample_id", "group"),
              [(sample, "ad" if sample in AD else "control") for sample in SAMPLES])

    baseline: dict[str, dict[str, float]] = {}
    for gene_number, gene in enumerate(GENES, 1):
        baseline[gene] = {}
        for sample_number, sample in enumerate(SAMPLES):
            value = 5.2 + (sample_number % 3) * 0.2 if gene_number % 7 == 0 else 7.0 + ((gene_number + sample_number) % 5) * 0.2
            baseline[gene][sample] = round(value, 2)
    write_tsv(raw / "expression.tsv", ("gene_id", "biotype", *SAMPLES),
              [(gene, "protein_coding", *(f"{baseline[gene][sample]:.2f}" for sample in SAMPLES)) for gene in GENES])

    splice_targets = [(sample, gene, f"Synthetic target {gene}") for index, sample in enumerate(AD)
                      for gene in (GENES[1], GENES[2], GENES[10 + index])]
    write_tsv(raw / "splicing_targets.tsv", ("sample_id", "gene_id", "gene_name"), splice_targets)
    ranks = {gene: {sample: 0.0 for sample in AD} for gene in GENES}
    for sample in AD:
        ranks[GENES[1]][sample] = 0.99  # Overlaps a splicing target.
        ranks[GENES[7]][sample] = 0.97
        ranks[GENES[8]][sample] = 0.96
        ranks[GENES[9]][sample] = 0.95  # Exact boundary must not be selected.
    write_tsv(raw / "pathogenicity_ranks.tsv", ("gene_id", "gene_symbol", *AD),
              [(gene, f"DEMO_{index:03d}", *(f"{ranks[gene][sample]:.2f}" for sample in AD))
               for index, gene in enumerate(GENES, 1)])
    write_tsv(raw / "reaction_subsystem.tsv", ("reaction_id", "subsystem"),
              [(reaction, PATHWAYS[index // 4]) for index, reaction in enumerate(REACTIONS)])

    states = [{gene: values.copy() for gene, values in baseline.items()}]
    audits: list[tuple[object, ...]] = []
    after_splicing = {gene: values.copy() for gene, values in states[-1].items()}
    for sample, gene, _ in splice_targets:
        before = after_splicing[gene][sample]
        after_splicing[gene][sample] = 0.0
        audits.append(("splicing", sample, gene, f"{before:.2f}", "0.00", "newly_zeroed"))
    states.append(after_splicing)
    after_both = {gene: values.copy() for gene, values in after_splicing.items()}
    for gene in GENES:
        for sample in AD:
            if ranks[gene][sample] > 0.95:
                before = after_both[gene][sample]
                after_both[gene][sample] = 0.0
                audits.append(("pathogenicity_rank_gt_0_95", sample, gene, f"{before:.2f}", "0.00",
                               "newly_zeroed" if before else "already_zero"))
    states.append(after_both)

    allowed = {(sample, gene) for sample, gene, _ in splice_targets}
    allowed.update((sample, gene) for gene in GENES for sample in AD if ranks[gene][sample] > 0.95)
    for condition_state in states[1:]:
        for gene in GENES:
            for sample in SAMPLES:
                if baseline[gene][sample] != condition_state[gene][sample] and (sample, gene) not in allowed:
                    raise RuntimeError(f"Unauthorised preview change: {sample}/{gene}")
                if sample in CONTROL and baseline[gene][sample] != condition_state[gene][sample]:
                    raise RuntimeError(f"Control value changed: {sample}/{gene}")
    if any(after_both[GENES[9]][sample] == 0 for sample in AD):
        raise RuntimeError("A rank exactly equal to 0.95 was incorrectly selected.")
    write_tsv(output / "audits" / "zeroing.tsv", ("layer", "sample_id", "gene_id", "value_before", "value_after", "status"), audits)

    activity: dict[str, dict[str, list[int]]] = {}
    summary_rows: list[tuple[object, ...]] = []
    sample_rows: list[tuple[object, ...]] = []
    for condition, values in zip(CONDITIONS, states):
        matrix = {sample: [int(values[gene][sample] >= THRESHOLD) for gene in GENES] for sample in SAMPLES}
        activity[condition] = matrix
        write_tsv(output / "results" / condition / "reaction_activity_matrix.tsv", ("reaction_id", *SAMPLES),
                  [(reaction, *(matrix[sample][index] for sample in SAMPLES)) for index, reaction in enumerate(REACTIONS)])
        for sample in SAMPLES:
            sample_rows.append((condition, sample, "ad" if sample in AD else "control", sum(matrix[sample])))
        for group, members in (("ad", AD), ("control", CONTROL)):
            counts = [sum(matrix[sample]) for sample in members]
            summary_rows.append((condition, group, len(members), f"{sum(counts) / len(counts):.2f}", min(counts), max(counts)))
    write_tsv(output / "results" / "sample_summary.tsv", ("condition", "sample_id", "group", "active_reactions"), sample_rows)
    write_tsv(output / "results" / "group_summary.tsv", ("condition", "group", "samples", "mean_active_reactions", "min", "max"), summary_rows)

    figures = output / "figures"
    bar_body = '<text x="36" y="42" font-size="24" font-weight="700">Active reaction counts</text><text x="36" y="67" font-size="13" fill="#60736f">Mean across four synthetic samples per group</text>'
    for index, (condition, label) in enumerate(zip(CONDITIONS, LABELS)):
        y = 110 + index * 115
        bar_body += f'<text x="36" y="{y}" font-size="15" font-weight="600">{html.escape(label)}</text>'
        for offset, group, color in ((20, "control", "#467d87"), (52, "ad", "#c26d55")):
            mean = float(next(row[3] for row in summary_rows if row[0] == condition and row[1] == group))
            length = mean * 21
            bar_body += f'<text x="36" y="{y + offset + 15}" font-size="13">{group.upper()}</text><rect x="116" y="{y + offset}" width="{length:.1f}" height="23" rx="5" fill="{color}"/><text x="{124 + length:.1f}" y="{y + offset + 17}" font-size="13">{mean:.2f} / 24</text>'
    svg(figures / "active_reactions.svg", 780, 465, "Active reaction counts by condition", "Synthetic thresholded reaction states for AD and control groups.", bar_body)

    heat_body = '<text x="34" y="42" font-size="23" font-weight="700">AD reaction activity by subsystem</text><text x="34" y="66" font-size="13">Fraction of synthetic AD samples with active reactions; descriptive preview</text>'
    for col, label in enumerate(LABELS):
        heat_body += f'<text x="{270 + col * 167}" y="112" font-size="12" text-anchor="middle">{html.escape(label)}</text>'
    for row, pathway in enumerate(PATHWAYS):
        y = 135 + row * 43
        heat_body += f'<text x="34" y="{y + 24}" font-size="13">{html.escape(pathway)}</text>'
        indices = range(row * 4, row * 4 + 4)
        for col, condition in enumerate(CONDITIONS):
            value = sum(activity[condition][sample][index] for sample in AD for index in indices) / 16
            shade = int(239 - value * 138)
            heat_body += f'<rect x="{210 + col * 167}" y="{y}" width="120" height="33" rx="5" fill="rgb({shade},{224 - int(value * 91)},{211 - int(value * 93)})"/><text x="{270 + col * 167}" y="{y + 22}" font-size="13" text-anchor="middle">{value:.0%}</text>'
    svg(figures / "subsystem_activity.svg", 780, 425, "AD reaction activity by subsystem", "Descriptive synthetic reaction activity, not pathway enrichment.", heat_body)

    target_body = '<text x="34" y="42" font-size="23" font-weight="700">Allowed expression zeroing</text><text x="34" y="67" font-size="13">Only synthetic AD targets were changed; controls were unchanged</text>'
    for index, sample in enumerate(AD):
        y = 112 + index * 69
        target_body += f'<text x="34" y="{y + 17}" font-size="14">{sample}</text><rect x="180" y="{y}" width="210" height="25" rx="4" fill="#467d87"/><rect x="390" y="{y}" width="140" height="25" rx="4" fill="#c26d55"/><text x="545" y="{y + 18}" font-size="13">3 splicing + 2 additional</text>'
    target_body += '<rect x="34" y="405" width="14" height="14" fill="#467d87"/><text x="56" y="417" font-size="13">Splicing targets</text><rect x="200" y="405" width="14" height="14" fill="#c26d55"/><text x="222" y="417" font-size="13">Additional pathogenicity targets</text>'
    svg(figures / "zeroing_audit.svg", 780, 445, "Synthetic expression zeroing audit", "Three splicing and two additional pathogenicity targets per AD sample.", target_body)

    rows = "".join(f'<tr><th scope="row">{html.escape(LABELS[CONDITIONS.index(condition)])}</th><td>{group.upper()}</td><td>{mean}</td></tr>'
                   for condition, group, _, mean, _, _ in summary_rows)
    report = f'''<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>Synthetic preview | Splicing and pathogenicity iMAT workflow</title>
<style>
:root{{--ink:#18332f;--muted:#526b66;--line:#cbdcd6;--bg:#f5f8f6;--accent:#236b70}}
*{{box-sizing:border-box}}body{{margin:0;background:var(--bg);color:var(--ink);font:16px/1.65 system-ui,-apple-system,Segoe UI,sans-serif}}
main{{max-width:1020px;margin:auto;padding:48px 22px 80px}}h1{{font-size:clamp(2rem,5vw,3.4rem);line-height:1.15;margin:.3rem 0 1rem}}h2{{margin-top:2.5rem}}
.eyebrow{{color:var(--accent);font-size:.8rem;letter-spacing:.12em;font-weight:700;text-transform:uppercase}}
.notice{{padding:16px 20px;border-left:5px solid #bc7654;background:#fff4e9;border-radius:0 8px 8px 0}}
.grid{{display:grid;grid-template-columns:repeat(auto-fit,minmax(260px,1fr));gap:18px}}figure{{margin:0;background:#fff;padding:15px;border:1px solid var(--line);border-radius:12px}}
img{{width:100%;height:auto}}figcaption{{font-size:.9rem;color:var(--muted)}}table{{width:100%;border-collapse:collapse;background:#fff}}
th,td{{border-bottom:1px solid var(--line);padding:10px;text-align:left}}th{{font-weight:650}}a{{color:var(--accent)}}a:focus-visible{{outline:3px solid #b76b42;outline-offset:3px}}
</style></head><body><main><p class="eyebrow">Open-source workflow preview</p><h1>Splicing and pathogenicity-informed reaction states</h1>
<p>This static report lets visitors inspect a small, reproducible input/output example before installing the scientific stack.</p>
<p class="notice"><strong>Synthetic preview only.</strong> These reaction states come from a simple threshold rule over invented values. This preview does <strong>not</strong> run Human-GEM, COBAMP, TROPPO iMAT, Fisher testing, or reporter-like statistics. It does not support biological or clinical conclusions.</p>
<h2>What the preview demonstrates</h2><p>The same four invented AD and four invented control samples pass through three conditions. Splicing targets are zeroed only in AD samples; genes with synthetic pathogenicity rank strictly greater than 0.95 are additionally zeroed. A target at exactly 0.95 is excluded. Every change is recorded in an <a href="audits/zeroing.tsv">audit table</a>.</p>
<div class="grid"><figure><a href="figures/active_reactions.svg" target="_blank" rel="noopener" aria-label="Open active reaction counts at full size"><img src="figures/active_reactions.svg" alt="Mean active reaction counts for synthetic AD and control groups across three conditions"></a><figcaption>Group means from toy reaction states. <a href="figures/active_reactions.svg" target="_blank" rel="noopener">View full size</a></figcaption></figure><figure><a href="figures/subsystem_activity.svg" target="_blank" rel="noopener" aria-label="Open subsystem activity at full size"><img src="figures/subsystem_activity.svg" alt="Heatmap of synthetic AD reaction activity fractions across six subsystem labels"></a><figcaption>Descriptive activity fractions, not enrichment. <a href="figures/subsystem_activity.svg" target="_blank" rel="noopener">View full size</a></figcaption></figure><figure><a href="figures/zeroing_audit.svg" target="_blank" rel="noopener" aria-label="Open expression zeroing audit at full size"><img src="figures/zeroing_audit.svg" alt="Three splicing and two additional pathogenicity targets zeroed per synthetic AD sample"></a><figcaption>Control expression remains unchanged. <a href="figures/zeroing_audit.svg" target="_blank" rel="noopener">View full size</a></figcaption></figure></div>
<h2>Generated summary</h2><table><thead><tr><th scope="col">Condition</th><th scope="col">Group</th><th scope="col">Mean active reactions / 24</th></tr></thead><tbody>{rows}</tbody></table>
<h2>Inspect the data</h2><ul><li><a href="data/raw/metadata.tsv">Sample metadata</a> and <a href="data/raw/expression.tsv">synthetic expression</a></li><li><a href="data/raw/splicing_targets.tsv">Splicing targets</a> and <a href="data/raw/pathogenicity_ranks.tsv">pathogenicity ranks</a></li><li><a href="results/group_summary.tsv">Group summary</a> and <a href="results/sample_summary.tsv">sample summary</a></li><li><a href="results/expression_only/reaction_activity_matrix.tsv">Expression-only matrix</a>, <a href="results/splicing_informed/reaction_activity_matrix.tsv">splicing matrix</a>, and <a href="results/splicing_pathogenicity_informed/reaction_activity_matrix.tsv">combined matrix</a></li></ul>
<p>Regenerate this report and every linked file with <code>python3 preview_demo.py</code>. The licensed full Human-GEM workflow is documented separately in the project README.</p>
</main></body></html>
'''
    (output / "index.html").write_text(report, encoding="utf-8")
    print(f"Wrote synthetic preview to {output}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT / "showcase", help="directory for generated preview files")
    args = parser.parse_args()
    generate(args.output.resolve())


if __name__ == "__main__":
    main()
