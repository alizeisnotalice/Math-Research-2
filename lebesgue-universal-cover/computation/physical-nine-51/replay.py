"""Verify the nine independent inputs at .83751 against the exterior and four native forests."""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
from replay_exterior import HERE, run, unpack


def equal(a, b):
    assert a.read_bytes() == b.read_bytes(), (a.name, b.name)


def forest(folder, driver, steps, model, row_file, expected):
    check = folder / "checked.json"
    arguments = [folder / "index.json", folder / "roots.json", folder / "stage0.jsonl",
                 check, "--initial-model", model]
    for i in range(1, steps + 1):
        stage = folder / f"stage{i}"
        arguments.extend(["--continuation", stage / "boxes.json", stage / row_file])
    run(driver, *arguments)
    equal(check, folder / expected)
    return check


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("computation", type=Path)
    parser.add_argument("tree", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--native-last", nargs=2, metavar=("TABLE64", "CLANG_CACHE"), type=Path,
                        help="Recompile the anchored kernel and visit all 128 final child boxes.")
    parser.add_argument("--workers", type=int, default=4)
    args = parser.parse_args()
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=False)
    comp, tree = args.computation.resolve(), args.tree.resolve()
    artifacts = json.loads((HERE / "result.json").read_bytes())
    for name, digest in artifacts["package_sha256"].items():
        assert sha256((HERE / name).read_bytes()).hexdigest() == digest, name
    ext = out / "exterior"
    run(HERE / "replay_exterior.py", comp, tree, ext)
    for name in ("core-2348", "core-1693", "core-1401", "geometry", "finite-788"):
        unpack(name + ".zip", out / name)
    g = out / "geometry"
    kit = g / "kit"
    base = out / "banks"
    points = [kit / f for f in ("p128.json", "physical-q.json", "distinct-q.json", "p15.json")]
    run(kit / "prepare.py", kit, base)
    run(kit / "verify_point_rows.py", kit / "nine-bank.json.gz", out / "coarse.json", *points)
    equal(out / "coarse.json", HERE / "coarse-angle-rows.json")
    fine = out / "fine"
    run(g / "refine_p128.py", base / "angular-bank.json", points[0], g / "slice/slice.cpp", fine,
        "--keep-nine")
    assert sha256((fine / "angular-nine-fine-p128.json").read_bytes()).hexdigest() == artifacts["fine_bank_sha256"]
    small = out / "core-1693"
    local = small / "local"
    combined = local / "combined220.jsonl"
    run(small / "merge_real_continuations.py", local / "stage2/boxes.json", combined,
        "--scope", local / "stage2a/boxes.json", local / "stage2a/rows.jsonl",
        "--scope", local / "stage2b/boxes.json", local / "stage2b/rows.jsonl")
    equal(combined, local / "stage2/rows.jsonl")
    north, south = out / "core-2348", out / "core-1401"
    forests = [
        forest(north, north / "audit_core51.py", 5, "seven", "result.jsonl", "final-audit.json"),
        forest(small / "prefix", small / "audit_core51.py", 4, "nine", "rows.jsonl", "final.json"),
        forest(local, small / "audit_core51.py", 2, "nine", "rows.jsonl", "final.json"),
        forest(south, south / "audit_core51.py", 4, "nine", "rows.jsonl", "final.json"),
    ]
    run(south / "map_last_seven51.py", south / "roots.json", south / "stage3/boxes.json",
        south / "stage3/rows.jsonl", south / "stage4/boxes.json", south / "stage4/rows.jsonl",
        forests[3], out / "last-seven.json")
    equal(out / "last-seven.json", south / "last-seven.json")
    old = comp / "global-083750/records/retained.json"
    theorem = out / "theorem.json"
    run(HERE / "join_51.py", old, ext / "exterior.json", *points, theorem, *forests)
    assert sha256(theorem.read_bytes()).hexdigest() == artifacts["theorem_sha256"]
    mapping = out / "map1024.json"
    map_args = []
    for k in ("a", "b", "c"):
        map_args.extend(["--replay", ext / ("layers-" + k)])
    for group in ("group01", "group2", "group3"):
        map_args.extend(["--direct-audit", ext / (group + ".json")])
    run(HERE / "map_1024_51.py", tree, ext / "seed/seed.json", ext / "layers/audit.json",
        old, theorem, mapping, *map_args)
    assert sha256(mapping.read_bytes()).hexdigest() == artifacts["canonical_map_sha256"]
    samples = out / "finite-788"
    run(samples / "finite51.py", samples / "inputs", out / "finite-check")
    expected = json.loads((samples / "finite/finite788.json").read_bytes())
    actual = json.loads((out / "finite-check/finite788.json").read_bytes())
    assert expected["groups"] == actual["groups"]
    assert (out / "finite-check/nine788.json.gz").read_bytes() == (samples / "finite/nine788.json.gz").read_bytes()
    run(HERE / "arc_thresholds.py", out / "arc-thresholds.json")
    equal(out / "arc-thresholds.json", HERE / "arc-thresholds.json")
    if args.native_last:
        table, cache = [i.resolve() for i in args.native_last]
        library = out / "anchor.so"
        import shutil
        subprocess.run([shutil.which("clang++") or "clang++", "-O3", "-std=c++17", "-fPIC", "-shared",
                        "-fsanitize=signed-integer-overflow", "-fsanitize-undefined-trap-on-error",
                        "-fmodules-cache-path=" + str(cache), "-I", str(g / "slice"), "-I", str(g / "include"),
                        str(g / "anchor-nine.cpp"), "-o", str(library)], check=True)
        index = out / "last-index.json"
        ids = sorted(int(x["leaf"].split("/")[1]) for x in json.loads((south / "stage4/boxes.json").read_bytes()))
        index.write_text(json.dumps({"1205": ids}) + "\n")
        result = out / "last-native.jsonl"
        # The driver imports its sibling sources before the repository's geometry routines.
        driver = g / "slice/run_slice.py"
        sys_path = f"import sys,runpy;sys.path[:0]=[{str(kit)!r},{str(driver.parent)!r}];sys.argv.pop(0);runpy.run_path(sys.argv[0],run_name='__main__')"
        subprocess.run([sys.executable, "-B", "-c", sys_path, str(driver), str(comp), str(tree),
                        str(table), str(library), str(fine / "angular-nine-fine-p128.json"), str(result),
                        "--target", ".83751", "--index", str(index), "--frontiers", str(south / "stage4/boxes.json"),
                        "--caps", "--mesh", "--cap", "4194303", "--workers", str(args.workers)], check=True)
        equal(result, south / "stage4/rows.jsonl")
    print("All 1,024 parts and 5,442 core leaves verified; no open boxes.")


if __name__ == "__main__":
    main()
