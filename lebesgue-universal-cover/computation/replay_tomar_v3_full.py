"""Replay all original roots from empty output and regenerate the public verification archive."""
import argparse
import json
from pathlib import Path
import shutil
import subprocess
import sys

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("checkout", type=Path, help="Tomar's v9 Git checkout")
parser.add_argument("certdir", type=Path, help="extracted frozen cert11d_08366 directory")
parser.add_argument("work", type=Path, help="new output directory")
parser.add_argument("--workers", default=2, type=int)
parser.add_argument("--sample-part", type=int, help="test this pipeline on one root without generating an all-roots claim")
args = parser.parse_args()
assert __debug__ and args.work.parent.is_dir() and args.work.parent.resolve().is_dir()
assert 1 <= args.workers <= 32
if args.sample_part is not None and not 0 <= args.sample_part < 2048:
    parser.error("sample part must lie between 0 and 2047")
args.work.mkdir(exist_ok=False)
work = args.work.resolve()
checkout, certdir = args.checkout.resolve(), args.certdir.resolve()
source = Path(__file__).resolve().parent

def run(program, *arguments, destination=None):
    command = [sys.executable, "-I", "-B", str(source/program), *map(str,arguments)]
    if destination is None:
        subprocess.run(command, check=True)
        return
    target = work/destination
    with target.with_suffix(".partial").open("w") as log:
        subprocess.run(command, stdout=log, check=True)
    target.with_suffix(".partial").replace(target)

table, library, native = work/"witnesses.npz", work/"integer.dylib", work/"integer"
scope = [] if args.sample_part is None else ["--parts", args.sample_part]
run("tomar_fan_table.py", certdir, table, destination="table.json")
run("run_tomar_integer_parts.py", certdir, native, library, table, "--fresh",
    "--workers", args.workers, *scope, destination="integer.json")
selected = range(2048) if args.sample_part is None else [args.sample_part]
assert {f.stem for f in native.glob("part*.json")} == {f"part{p:04}" for p in selected}
final = work/"arb"
final.mkdir()
todo = list(selected)
for depth in (12,24,32,40,48):
    current = work/f"arb-depth-{depth}"
    restricted = [] if len(todo)==2048 else ["--parts", ",".join(map(str,todo))]
    run("finish_tomar_v3_parts.py",checkout,certdir,native,current,"--depth",depth,*restricted,
        destination=f"arb-depth-{depth}.json")
    pending = []
    for part in todo:
        name = f"part{part:04}.json"
        receipt = json.loads((current/name).read_bytes())
        if receipt["unresolved"]:
            pending.append(part)
        else:
            shutil.copyfile(current/name,final/name)
    todo = pending
    if not todo:
        break
if todo:
    raise RuntimeError(f"roots still missing area checks: {todo}")
pair = ["--pair","local",native,final]
strict = ["--require-complete"] if args.sample_part is None else []
run("audit_tomar_v3_coverage.py",checkout,certdir,table,*pair,*strict,destination="coverage.json")
ledger = json.loads((work/"coverage.json").read_bytes())
assert ledger["certified_global_11d_replay"] == (args.sample_part is None)
run("bundle_tomar_v3_receipts.py",work/"coverage.json",work/"private-receipts.zip",*pair,
    destination="archive.json")
run("publish_tomar_v3_archive.py",work/"private-receipts.zip",work/"receipts.zip",destination="publication.json")
public_scope=["--partial"] if args.sample_part is not None else []
run("audit_tomar_v3_public.py",checkout,certdir,table,work/"receipts.zip",work,*public_scope,destination="public-audit.json")
if args.sample_part is None:
    run("summarize_tomar_v3_archive.py",work/"receipts.zip",destination="summary.json")
print(json.dumps({k:ledger[k] for k in ("target","certified_global_11d_replay","certified_root_count",
                                       "original_leaves","integer_closed","arb_closed","arb_failed")},indent=2))
