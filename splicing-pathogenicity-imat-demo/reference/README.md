# Human-GEM 1.12.0 reference

The supplied `Human-GEM_1.12.0.xml` has SHA-256
`53a74c4495ca4af1334ce67b74cc09fc1ab30c7648eef561f478cac2363f13ab`.
It is byte-identical to the pinned
[Human-GEM v1.12.0 release model](https://github.com/SysBioChalmers/Human-GEM/releases/tag/v1.12.0).
The full XML and upstream release tree remain local because the public
repository supplies `download_model.py` to retrieve and verify the exact model.

`model_index.tsv.gz` is generated from that verified XML by
`build_model_index.py`. It contains 13,070 reaction IDs, their Human-GEM
subsystem memberships, and linked model gene IDs. One reaction has two
subsystem memberships, yielding 13,071 reaction–subsystem rows when expanded.
The index is reference annotation only: it is **not** an iMAT output or a
reaction-activity matrix.

Human-GEM is distributed under
[CC BY 4.0](../third_party_licenses/Human-GEM-CC-BY-4.0.md). Cite the upstream
project when using the model or its derived annotations.
