# Human-GEM model source

The full pipeline uses `Human-GEM-1.12.0.xml`, the SBML model from the
[Human-GEM v1.12.0 release](https://github.com/SysBioChalmers/Human-GEM/releases/tag/v1.12.0).
The model binary is excluded from this repository to keep the project small.

Download and verify it with:

```bash
python3 download_model.py
```

The downloader uses the version-pinned upstream file and checks SHA-256:
`53a74c4495ca4af1334ce67b74cc09fc1ab30c7648eef561f478cac2363f13ab`.
Human-GEM is licensed under CC BY 4.0; see the copy in
`third_party_licenses/Human-GEM-CC-BY-4.0.md` and cite the original project.
