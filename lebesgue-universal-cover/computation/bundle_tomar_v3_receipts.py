"""Archive exactly the audited receipts in a V3 coverage ledger with stable ZIP names."""
import argparse
import hashlib
import json
from pathlib import Path
from zipfile import ZIP_DEFLATED, ZipFile, ZipInfo

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("ledger", type=Path)
parser.add_argument("output", type=Path)
parser.add_argument("--pair", nargs=3, action="append", required=True, metavar=("HOST","NATIVE","ARB"))
args = parser.parse_args()
ledger = json.loads(args.ledger.read_bytes())
roots = {}
for host, native, arb in args.pair:
    for file in Path(native).glob("part*.json"):
        key = (host, file.stem)
        assert key not in roots
        roots[key] = (file, Path(arb) / file.name)
data = []
expected = {}
for row in ledger["rows"]:
    part, host = row[:2]
    stem = f"part{part:04}"
    native, supplement = roots[(host, stem)]
    for name, file, sha in zip(("integer", "integer", "arb"), (native,native.with_suffix(".npz"),supplement), row[-3:]):
        label = f"{host}/{name}/{file.name}"
        assert label not in expected and hashlib.sha256(file.read_bytes()).hexdigest() == sha
        expected[label] = sha
        data.append((label, file))
data.append(("coverage.json", args.ledger))
expected["coverage.json"] = hashlib.sha256(args.ledger.read_bytes()).hexdigest()
tmp = args.output.with_suffix(".partial")
with ZipFile(tmp,"w",compression=ZIP_DEFLATED,compresslevel=9) as archive:
    for label, file in data:
        entry = ZipInfo(label)
        entry.external_attr = 0o100644 << 16
        archive.writestr(entry,file.read_bytes(),compress_type=ZIP_DEFLATED,compresslevel=9)
with ZipFile(tmp) as archive:
    assert set(archive.namelist()) == set(expected) and len(archive.namelist()) == len(expected)
    for label, sha in expected.items():
        assert hashlib.sha256(archive.read(label)).hexdigest() == sha
tmp.replace(args.output)
print(json.dumps({"parts":len(ledger["rows"]), "objects":len(expected), "bytes":args.output.stat().st_size,
                  "sha256":hashlib.sha256(args.output.read_bytes()).hexdigest()},indent=2))
