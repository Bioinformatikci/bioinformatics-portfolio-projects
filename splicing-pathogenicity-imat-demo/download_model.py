"""Download the exact Human-GEM 1.12.0 SBML model used by the full pipeline."""
from __future__ import annotations

import argparse
import hashlib
import shutil
import tempfile
import urllib.request
from pathlib import Path


ROOT = Path(__file__).resolve().parent
TARGET = ROOT / "model" / "Human-GEM-1.12.0.xml"
SOURCE = "https://raw.githubusercontent.com/SysBioChalmers/Human-GEM/v1.12.0/model/Human-GEM.xml"
SHA256 = "53a74c4495ca4af1334ce67b74cc09fc1ab30c7648eef561f478cac2363f13ab"


def digest(path: Path) -> str:
    hasher = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            hasher.update(chunk)
    return hasher.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--overwrite", action="store_true", help="replace an existing mismatched model file")
    args = parser.parse_args()
    if TARGET.exists():
        current = digest(TARGET)
        if current == SHA256:
            print(f"Model already verified: {TARGET}")
            return
        if not args.overwrite:
            raise SystemExit(f"Existing model has SHA-256 {current}; pass --overwrite to replace it.")
    TARGET.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(prefix="human-gem-", suffix=".xml", dir=TARGET.parent, delete=False) as stream:
        temporary = Path(stream.name)
        try:
            request = urllib.request.Request(SOURCE, headers={"User-Agent": "splicing-pathogenicity-imat-demo/1.0"})
            with urllib.request.urlopen(request, timeout=60) as response:
                shutil.copyfileobj(response, stream)
        except BaseException:
            temporary.unlink(missing_ok=True)
            raise
    actual = digest(temporary)
    if actual != SHA256:
        temporary.unlink(missing_ok=True)
        raise SystemExit(f"Downloaded model checksum mismatch: {actual}")
    temporary.replace(TARGET)
    print(f"Downloaded and verified Human-GEM v1.12.0: {TARGET}")


if __name__ == "__main__":
    main()
