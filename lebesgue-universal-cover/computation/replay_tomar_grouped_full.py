"""Run the native and deep D3 five-body checks from empty output, then audit and seal all results."""
import argparse
from pathlib import Path
import subprocess
import sys

p = argparse.ArgumentParser(description=__doc__)
for key in ("certdir", "slice_index", "work", "cache"):
    p.add_argument(key, type=Path)
p.add_argument("--workers", type=int, default=1)
p.add_argument("--sample-parts", help="comma-separated slice roots for a partial pipeline test")
args = p.parse_args()
assert __debug__ and 1 <= args.workers <= 32
args.work.mkdir(exist_ok=False)
work = args.work.resolve()
here = Path(__file__).resolve().parent

def run(program, *argv, name):
    path = work / name
    with path.with_suffix(".partial").open("w") as out:
        subprocess.run([sys.executable, "-I", "-B", here/program, *map(str, argv)], check=True, stdout=out)
    path.with_suffix(".partial").replace(path)

table = work / "witnesses32.npz"
native, deep = work/"native", work/"deep"
scope = [] if args.sample_parts is None else ["--parts", args.sample_parts]
partial = [] if args.sample_parts is None else ["--partial"]
run("tomar_fan_table.py", args.certdir, table, "--subdiv", 32, name="table.json")
run("run_tomar_grouped_slice.py", args.certdir, native, work/"libtomar_grouped_scored.dylib", table,
    "--workers", args.workers, "--slice-index", args.slice_index, "--cache", args.cache, *scope, name="native.json")
run("check_tomar_grouped_residues.py", args.certdir, table, native, deep, args.cache, name="deep.log")
ledger = work / "audit.json"
run("audit_tomar_grouped_slice.py", args.certdir, table, args.slice_index, native, deep, work, args.cache, *partial,
    name=ledger.name)
public = work/"public"
run("bundle_tomar_grouped_slice.py", ledger, native, deep, public, name="publication.json")
run("audit_tomar_grouped_public.py", args.certdir, table, args.slice_index, ledger, public/"tomar-grouped-08371-index.json",
    work, args.cache, *partial, name="public-audit.json")
print((work/"public-audit.json").read_text())
