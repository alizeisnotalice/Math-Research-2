"""Replay the historical six-set forest and the scaled nine-set lower bound."""
import argparse
from fractions import Fraction
import gzip
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import zipfile

HERE = Path(__file__).resolve().parent
C = HERE.parent
DIGEST = {
    "six": "a757bf5c172bf1909d88c3b65156a6e958a629866c9a3ef7b384fadac52c8a50",
    "scale": "c2b030544c6565a8379aaa271cec35582d060486463866f8dabeac7ccf7b90ea",
    "closest": "6e893ff6ce9236858112187d9a97b6dc299dd5bd6913f2db140f701121e9587d",
    "fine-501": "e85192f362d312475bede0351886c39b44a06cca769267af30e669f75c45a8aa",
    "independent-scale": "80e8b7cbc64cd8dfbc1d99fafa24bd438ec5c84d2ffb058289e6f65c501955e8",
}


def unpack(archive, folder):
    with zipfile.ZipFile(archive) as bundle:
        for name in bundle.namelist():
            path = Path(name)
            assert not path.is_absolute() and ".." not in path.parts
        bundle.extractall(folder)


def read_json(path):
    return json.loads(path.read_bytes())


def run(script, *params):
    subprocess.run([sys.executable, "-B", str(script), *map(str, params)],
                   cwd=script.parent, check=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("tree", "table64", "cache", "merged45", "work"):
        parser.add_argument(name, type=Path)
    parser.add_argument("--native-501", action="store_true", help="rerun the 501 primary integer trees")
    parser.add_argument("--six-native", choices=["last", "all"], help="rebuild the six-body integer trees")
    parser.add_argument("--workers", type=int, default=3)
    args = parser.parse_args()
    assert __debug__ and args.workers > 0
    work = args.work.resolve()
    tree, table, cache, merged = (x.resolve() for x in (args.tree, args.table64, args.cache, args.merged45))
    work.mkdir(parents=True, exist_ok=False)
    for name, expected in DIGEST.items():
        archive = HERE / f"{name}.zip"
        assert hashlib.sha256(archive.read_bytes()).hexdigest() == expected
        unpack(archive, work / name)
    six = work / "six"
    bank = six / "fine30.json"
    bank.write_bytes(gzip.decompress((six / "fine30.json.gz").read_bytes()))
    run(six / "verify_original_bank.py", C, bank, six / "fresh-support.json")
    assert read_json(six / "fresh-support.json") == read_json(six / "physical-support-audit.json")
    unpack(C / "whole-leaves-083750-second/wave2.zip", work / "r248")
    forest = six / "forest.json.gz"
    forest_args = [C, tree, work / "r248/public-replay-r2-083750",
                   six / "ledger-eight.json", six / "selectionOld30.json.gz", forest,
                   six / "six-five-floors.json", "--merged45", merged]
    if args.six_native:
        build = work / "six-build"
        run(six / "build_six_kernels.py", forest, build, "--extra=-fmodules-cache-path=" + str(cache))
        forest_args += ["--native", "--src", six, "--table", table, "--bank", bank,
                        "--libdir", build, "--scratch", work / "six-kernel", "--workers", args.workers]
        if args.six_native == "last":
            forest_args += ["--only", "1205/152529"]
    run(six / "audit_forest.py", *forest_args)
    six_data = json.loads(gzip.decompress(forest.read_bytes()))
    assert not six_data["open"] and six_data["baseline"] + len(six_data["roots"]) == 5442
    six_lower = Fraction(six_data["integer_goal"], 2**57)
    assert six_lower - Fraction(".83750") == Fraction(1, 720575940379279360)

    primary = work / "independent-scale/frozen/p66-oldq-501.jsonl"
    scale = work / "scale"
    out = work / "normalized"
    args_scale = [C / "seven-compression-61/scaled.zip", C, tree, table, cache, out,
                  "--workers", args.workers]
    if not args.native_501:
        args_scale += ["--reuse-roots", primary]
    run(scale / "replay.py", *args_scale)
    if args.native_501:
        assert primary.read_bytes() == (out / "roots51.jsonl").read_bytes()
    fine = work / "fine-501"
    run(fine / "audit_501.py", fine / "native/keys501.json", primary,
        fine / "native/fine-nine-501.jsonl", work / "two-501.json", "--require-fine-complete")
    assert read_json(work / "two-501.json") == read_json(fine / "audit-501.json")
    for n in (730, 1931, 1991):
        actual = out / f"finite{n}/nine-finite.json.gz"
        assert gzip.decompress(actual.read_bytes()) == gzip.decompress((HERE / f"finite{n}.json.gz").read_bytes())
    run(work / "closest/replay.py", work / "pair-replay")
    geom = read_json(out / "geometry51.json")
    nine_lower = Fraction(geom["G"]) / Fraction(geom["finite_scale"])**2
    assert nine_lower == Fraction(geom["source_area"])
    assert nine_lower > Fraction(".83751146695094269361337962811932078")
    result = {
        "historical_six_lower": str(six_lower), "historical_six_unresolved": len(six_data["open"]),
        "historical_six_previous": six_data["baseline"], "historical_six_added": len(six_data["roots"]),
        "scaled_nine_lower": str(nine_lower), "scaled_nine_decimal": geom["source_decimal"],
        "exceptional_roots": read_json(work / "two-501.json"),
        "individual_group_bounds": [read_json(out / f"finite{n}/finite-report.json")["finite_decimal"]
                                    for n in (730, 1931, 1991)],
        "closest_eligible_input_pair": ["oldQ", "newQ"],
        "pair_diameter_interval": ["1.00067004061890625", "1.000673064"],
    }
    (work / "result.json").write_text(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    main()
