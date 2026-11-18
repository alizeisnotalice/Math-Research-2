"""Check the shard coverage and verify both polynomial gcd kernels."""

import argparse
from hashlib import sha256
import json
from pathlib import Path
import random
import subprocess

import sympy as sp


def cover(directory):
    summary = json.loads((directory / "summary.json").read_bytes())
    assert summary["status"] == "complete" and summary["masks"] == 1 << 23
    assert summary["shards"] == 128 and summary["shard_size"] == 1 << 16
    numbered = sorted(directory.glob("[0-9]*.json"))
    assert len(numbered) == 128
    rows = []
    for i, part in enumerate(numbered):
        start = i << 16
        assert part.name == f"{start:08d}.json"
        row = json.loads(part.read_bytes())
        assert (row["start"], row["stop"]) == (start, start + (1 << 16))
        assert row["first_prime_pass"] == 0 and not row["hits"]
        row.pop("pid")
        rows.append(row)
    assert summary["first_prime_pass"] == 0 and not summary["hits"]
    return dict(summary=summary, shards=rows)


def fixture(work, source, name, code):
    fixture_source = work / f"{name}-fixture.cpp"
    fixture_source.write_text(
        f'#define main {name}_entry\n#include "{source / (name + ".cpp")}"\n'
        f'#undef main\nint main(){{{code};return 0;}}\n'
    )
    binary = fixture_source.with_suffix("")
    subprocess.run(["/usr/bin/clang++", "-O3", "-std=c++20", str(fixture_source), "-o", str(binary)], check=True)
    subprocess.run([str(binary)], check=True)


def main():
    assert __debug__
    parser = argparse.ArgumentParser(description=__doc__)
    for option in ("ordinary", "reciprocal", "bin_dir", "scratch", "out"):
        parser.add_argument(option, type=Path)
    args = parser.parse_args()
    first, second = cover(args.ordinary), cover(args.reciprocal)
    ordinary = args.bin_dir.resolve() / "scan64"
    reciprocal = args.bin_dir.resolve() / "scan64_reciprocal"
    source = Path(__file__).resolve().parent
    args.scratch.mkdir(exist_ok=False, parents=True)
    fixture(args.scratch, source, "scan64",
            "Polynomial a{};a[0]=a[64]=1;a[32]=161;"
            "if(gcd_degree<163>(a,true)!=32)return 2;"
            "a[32]=1277;if(gcd_degree<1279>(a,false)!=32)return 3")
    fixture(args.scratch, source, "scan64_reciprocal",
            "Poly a{};a[0]=a[64]=1;a[32]=-2;"
            "if(exact_degree<163>(a)!=32||exact_degree<1279>(a)!=32)return 2")
    seed = 20261118
    masks = [0, 1, 1 << 22, (1 << 23) - 1] + random.Random(seed).sample(range(2, (1 << 22) - 1), 64)
    z = sp.Symbol("x")
    fixed = sum(z**i for i in range(36) if i != 7) + z**38 + z**64
    powers = tuple(i for i in range(36, 60) if i != 38)
    for mask in masks:
        coeffs = fixed + sum(z**i for j, i in enumerate(powers) if mask & (1 << j))
        expected = []
        for p in (163, 1279):
            poly = sp.Poly(coeffs, z, modulus=p)
            expected.append(poly.gcd(poly.diff()).degree())
        for binary in (ordinary, reciprocal):
            assert json.loads(subprocess.check_output([str(binary), "--inspect", str(mask)])) == expected
    digests = {path.name: sha256(path.read_bytes()).hexdigest() for path in sorted(source.glob("scan64*.cpp"))}
    report = dict(status="verified", fixtures=4, sampled_masks=len(masks), seed=seed,
                  sympy=sp.__version__, source_sha256=digests, ordinary=first, reciprocal=second)
    with args.out.open("x") as out:
        json.dump(report, out, indent=2)
        out.write("\n")
    print(report["status"], len(masks))


if __name__ == "__main__":
    main()
