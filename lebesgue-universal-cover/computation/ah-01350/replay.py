"""Reassemble the frozen .01350 proof by file content and run its original verifier."""
import argparse
from hashlib import sha256
import lzma
from pathlib import Path
import subprocess
import sys
import tarfile

HERE = Path(__file__).resolve().parent
SOURCE_SHA = "69bcff542ce683117b1b534c1d66e734769cf753c382b539ebfe8fccc58f15ad"
RAW_ROWS = {
    "j32-34.json.xz": (
        "segments/j32_34/verified-32-34.json",
        "ce52e9254e1e961d1a672201da705fece33310e470a9a46ca95884114eae1861",
    ),
    "j35-38.json.xz": (
        "segments/j35_38/full-35-38.json",
        "9302e686810623f85221889b8e62b568d85001c13ffe91a63b10c11cd25db244",
    ),
}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("work", type=Path)
    parser.add_argument("--full", action="store_true", help="repeat every original interval kernel")
    args = parser.parse_args()
    work = args.work.resolve()
    work.mkdir(parents=True, exist_ok=False)
    kit = work / "kit"
    kit.mkdir()
    core = HERE / "source.tar.xz"
    if sha256(core.read_bytes()).hexdigest() != SOURCE_SHA:
        raise ValueError("Unexpected certificate archive.")
    with tarfile.open(core, "r:xz") as archive:
        members = archive.getmembers()
        if len(members) != 1005 or not all(m.isfile() for m in members):
            raise ValueError("Unexpected certificate members.")
        archive.extractall(kit, members=members, filter="data")
    for file_name, (relative, digest) in RAW_ROWS.items():
        raw = lzma.decompress((HERE / file_name).read_bytes())
        if sha256(raw).hexdigest() != digest:
            raise ValueError(f"Unexpected interval rows for {relative}")
        (kit / relative).write_bytes(raw)
    output = work / "checked"
    command = [sys.executable, "-I", "-B", str(kit / "verify01350.py"), str(output)]
    if args.full:
        command.append("--full")
    subprocess.run(command, cwd=kit, check=True)
    actual = (output / "cover-result.json").read_bytes()
    if not args.full and actual != (HERE / "cover.json").read_bytes():
        raise ValueError("Replayed certificate differs from its frozen receipt.")


if __name__ == "__main__":
    main()
