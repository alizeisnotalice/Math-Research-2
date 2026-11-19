"""Replay the complete 0.01400 universal cover certificate from the public parts."""
import argparse
import hashlib
from pathlib import Path
import subprocess
import sys
import zipfile


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("work", type=Path, help="a new directory for extracted files and results")
    modes = parser.add_mutually_exclusive_group()
    modes.add_argument("--full", action="store_true", help="rebuild each original interval certificate")
    modes.add_argument("--assemble-only", action="store_true", help="only verify and reassemble the archives")
    args = parser.parse_args()
    package = Path(__file__).resolve().parent
    work = args.work.resolve()
    work.mkdir(parents=True, exist_ok=False)
    source = work / "source"
    archive = package / "source.zip"
    expected = "833328bdd10af144a46cd480d143ba162db8f30d330b744a8dcfaee7fdb276f6"
    if hashlib.sha256(archive.read_bytes()).hexdigest() != expected:
        raise ValueError("source.zip has a different SHA-256 digest")
    with zipfile.ZipFile(archive) as packed:
        packed.extractall(source)
    parts = [package / name for name in ("j32-35-a.zip", "j32-35-b.zip", "j32-35-c.zip", "j36-39.zip")]
    cmd = [sys.executable, "-I", "-B", source / "bootstrap014.py", *parts, work / "verification"]
    if args.full:
        cmd.append("--full")
    if args.assemble_only:
        cmd.append("--assemble-only")
    subprocess.run(cmd, check=True)
    if not (args.full or args.assemble_only):
        receipt = work / "verification" / "checked" / "result.json"
        if receipt.read_bytes() != (package / "result.json").read_bytes():
            raise ValueError("the recomputed receipt differs from the published receipt")


if __name__ == "__main__":
    main()
