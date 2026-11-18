"""Verify saved replay receipts against the original tree and frozen first-pass proofs."""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import re
import subprocess
import sys
import zipfile


def main(tree, work):
    package = Path(__file__).resolve().parent
    comp = package.parent.parent
    work.mkdir(parents=True, exist_ok=False)
    ordinary = []
    dense = 0
    with zipfile.ZipFile(package / "ordinary-replay.zip") as archive:
        names = archive.namelist()
        assert len(names) == len(set(names)) == 3095
        for name in names:
            match = re.fullmatch(r"(\d+)-(\d+)\.json", name)
            if not match:
                assert re.fullmatch(r"dense/\d+-\d+/result\.json", name)
                dense += 1
            path = work / name
            path.parent.mkdir(parents=True, exist_ok=True)
            data = archive.read(name)
            path.write_bytes(data)
            if match:
                ordinary.append([*map(int, match.groups()), str(path.resolve()), sha256(data).hexdigest()])
    assert (len(ordinary), dense) == (3071, 24)
    manifest = work / "ordinary-index.json"
    manifest.write_text(json.dumps(ordinary))
    ledger = work / "verified.json"
    subprocess.run([
        sys.executable, "-B", str(package / "full_leaf_ledger.py"), str(package.parent / "volumes"),
        str(comp), str(tree), str(work / "dense"), str(manifest), str(ledger),
        "--pure-replay", str(package.parent / "records" / "pure45.json"),
    ], check=True)
    assert ledger.read_bytes() == (package / "full-leaf-ledger.json").read_bytes()
    print("Verified 3273 terminal leaves and 3071 byte-identical ordinary receipts.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("tree", type=Path)
    parser.add_argument("work", type=Path)
    args = parser.parse_args()
    main(args.tree, args.work)
