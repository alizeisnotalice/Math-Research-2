#!/usr/bin/env python3
"""Portable launcher for the frozen round-046 verification programs.

The frozen programs retain their original bytes and therefore their original
relative-neighbour assumptions.  This launcher materializes those neighbours
under a caller-chosen temporary directory, never modifying the package.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def validate() -> None:
    manifest = json.loads((ROOT / "SHA256SUMS.json").read_text())
    failures = []
    for row in manifest["files"]:
        p = ROOT / row["path"]
        if not p.is_file() or p.stat().st_size != row["bytes"] or sha(p) != row["sha256"]:
            failures.append(row["path"])
    if failures:
        raise SystemExit("hash validation failed: " + ", ".join(failures))
    print(f"OK: {len(manifest['files'])} package files match SHA256SUMS.json")


def materialize(work: Path) -> Path:
    """Create the unmodified historical relative layout below *work*."""
    if work.exists() and any(work.iterdir()):
        raise SystemExit(f"work directory must be absent or empty: {work}")
    parent = work / "round046-parent"
    linked = work / "agents/coprime_stretches/round046/core6_linked"
    single = work / "agents/algebraic_global/round046/core6_single"
    parent.mkdir(parents=True, exist_ok=True)
    linked.parent.mkdir(parents=True, exist_ok=True)
    single.parent.mkdir(parents=True, exist_ok=True)
    for p in (ROOT / "frozen/parent").iterdir():
        shutil.copy2(p, parent / p.name)
    for name in ("linked-prefix-input.json", "linked-input-tuples.json"):
        shutil.copy2(ROOT / "data/parent" / name, parent / name)
    # The linked and single paths are read-only source/data trees; a symlink
    # keeps a replay from overwriting an archived result.
    linked.symlink_to(ROOT / "data/linked", target_is_directory=True)
    single.symlink_to(ROOT / "data/single", target_is_directory=True)
    return parent


def run(args: argparse.Namespace) -> None:
    parent = materialize(args.workdir)
    if args.mode == "linked-b":
        cmd = [sys.executable, str(parent / "verify_linked_full.py"), "--b", str(args.b)]
    elif args.mode == "single":
        # This also writes a fresh large terminal journal under --workdir.
        cmd = [sys.executable, str(parent / "verify_single_core.py")]
    elif args.mode == "direct-prefix":
        cmd = [sys.executable, str(parent / "direct_linked_prefix.py")]
    else:
        raise AssertionError(args.mode)
    print("running:", " ".join(cmd))
    subprocess.run(cmd, check=True)


def compare(args: argparse.Namespace) -> None:
    script = ROOT / "frozen/parent/compare_linked_trees.py"
    direct = ROOT / "data/parent/linked-direct-trees.jsonl.gz"
    trees = [ROOT / "data/linked/prefix_tree_0_3708.jsonl.gz", ROOT / "data/linked/prefix_tree_3708_7417.jsonl.gz"]
    cmd = [sys.executable, str(script), str(direct), *(str(x) for x in trees), "--out", str(args.out)]
    print("running:", " ".join(cmd))
    subprocess.run(cmd, check=True)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest="command", required=True)
    sub.add_parser("validate", help="check all archived package-file hashes")
    pc = sub.add_parser("compare-trees", help="recompare all archived direct/HGQ tree levels")
    pc.add_argument("--out", type=Path, required=True, help="output JSON outside the package")
    pr = sub.add_parser("replay", help="run one expensive frozen verifier in a fresh work directory")
    pr.add_argument("mode", choices=("linked-b", "single", "direct-prefix"))
    pr.add_argument("--workdir", type=Path, required=True, help="empty directory outside this package")
    pr.add_argument("--b", type=int, help="required only for linked-b; choose a nonzero b in [-7,7]")
    args = p.parse_args()
    if args.command == "validate":
        validate()
    elif args.command == "compare-trees":
        compare(args)
    else:
        if args.mode == "linked-b" and (args.b is None or args.b == 0 or abs(args.b) > 7):
            raise SystemExit("linked-b requires --b in {-7,...,-1,1,...,7}")
        if args.mode != "linked-b" and args.b is not None:
            raise SystemExit("--b is only valid for linked-b")
        run(args)


if __name__ == "__main__":
    main()
