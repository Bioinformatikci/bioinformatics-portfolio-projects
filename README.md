# Bioinformatics Portfolio Projects

This repository collects selected bioinformatics portfolio projects, analysis scripts, workflow notes, and small applications from my computational biology training.

The projects focus on sequence analysis, transcriptomics, alternative splicing workflows, exploratory visualization, and compact reproducible analysis utilities.

## Repository Overview

| Folder | Focus | Main tools |
| --- | --- | --- |
| `sequence-alignment-and-scoring-matrices/` | PAM scoring matrices and sequence alignment methods | Python, NumPy, pandas, matplotlib, Jupyter |
| `protein-domain-annotation-mini-organism/` | ORF-level protein annotation with BLAST, Pfam, and marker-gene placement | Python, Biopython, pandas, matplotlib |
| `geo-psoriasis-expression-analysis/` | GEO psoriasis expression analysis with differential testing and result plots | R, GEOquery, Biobase |
| `bioinformatics-python-scripts/` | Small reproducible bioinformatics utilities and analysis scripts | Python, pandas, SciPy, matplotlib |
| `rna-seq-alternative-splicing-rmats/` | RNA-seq alternative splicing workflow notes | STAR, rMATS, rmats2sashimiplot, shell |
| `interactive-shiny-dashboards/` | Interactive data exploration applications | R, Shiny, ggplot2, DT |
| [`splicing-pathogenicity-imat-demo/`](splicing-pathogenicity-imat-demo/) | Three-stage Human-GEM iMAT workflow on a synthetic 200-sample cohort, with computed matrices, figures, and an HTML results report | Python, COBAMP, TROPPO, Human-GEM, Gurobi, R |

## Highlighted Work

- Sequence alignment notebook with PAM-style matrix calculations and a rendered HTML version for quick review.
- Protein-domain annotation project with ORF annotations, Pfam architecture plots, and marker-gene phylogenetic placement.
- GEO psoriasis expression analysis with reproducible R code, summary tables, and diagnostic plots.
- Python scripts for FASTA parsing, motif search, simple genome assembly logic, GBM expression summaries, volcano plotting, and variant population summaries.
- R/Shiny applications for interactive Palmer Penguins and obesity transcriptome data exploration.
- RNA-seq alternative splicing workflow notes using STAR, rMATS, and rmats2sashimiplot.
- [Synthetic splicing and pathogenicity iMAT workflow](splicing-pathogenicity-imat-demo/): inspect the [vector figures](splicing-pathogenicity-imat-demo/showcase/figures/active_reactions.svg), [computed reaction matrices and statistics](splicing-pathogenicity-imat-demo/showcase/results/), [supplied synthetic inputs](splicing-pathogenicity-imat-demo/synthetic_single_cohort/), and [downloadable HTML report](splicing-pathogenicity-imat-demo/showcase/index.html). The project README documents the verified Human-GEM release, analysis steps, solver requirements, and reproduction limits. Synthetic differences are not biological or clinical findings.

![Computed iMAT reaction activity in a synthetic cohort](splicing-pathogenicity-imat-demo/showcase/figures/active_reactions.svg)

## Topics Covered

- Sequence alignment and substitution scoring matrices
- Protein domain annotation and homology-based functional interpretation
- Public GEO expression-set analysis
- FASTA parsing and motif search
- Differential expression-style summary analysis
- Variant frequency summaries and visualization
- Alternative splicing workflow design with rMATS
- Interactive R/Shiny visualization

## Notes On Authorship

The repository contains Burak Keskin's portfolio work. The synthetic iMAT demonstration includes analysis code, artificial inputs, and presentation materials developed with Codex assistance. Its README identifies the calculated outputs, source model, reproducibility requirements, and interpretation limits.

## Repository Status

This is a portfolio-oriented repository. Some scripts are compact analysis examples or workflow notes rather than fully packaged software. Where external datasets are required, the relevant project README explains the expected input files.

Large raw sequencing files, model archives, and temporary analysis outputs are excluded from version control. The synthetic iMAT project includes curated compressed reaction matrices and result tables so the figures can be audited without rerunning the solver. Small coursework datasets needed by the Shiny examples are included so those apps can be reviewed locally.

## Suggested Use

Start with the project-level README files inside each folder. They explain what each project contains, which tools are involved, and whether the files are intended to be run directly or reviewed as coursework/workflow examples.
