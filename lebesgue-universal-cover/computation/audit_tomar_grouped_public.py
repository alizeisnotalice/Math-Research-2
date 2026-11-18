"""Unpack all hash-bound public grouped-fan archives against the authenticated root ledger."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
from zipfile import ZipFile

def sha(data):
    return hashlib.sha256(data).hexdigest()

def main():
    p = argparse.ArgumentParser(description=__doc__)
    for key in ("certdir", "table", "slice_index", "ledger", "bundle_index", "scratch", "cache"):
        p.add_argument(key, type=Path)
    p.add_argument("--partial", action="store_true")
    a = p.parse_args()
    index = json.loads(a.bundle_index.read_bytes())
    ledger = json.loads(a.ledger.read_bytes())
    assert sha(a.ledger.read_bytes()) == index["audit_sha256"]
    assert index["target"] == ledger["target"] == "0.8371"
    assert index["parts"] == ledger["parts"] and index["complete"] == ledger["complete"] == (not a.partial)
    expected = {r[0]: r[1:] for r in ledger["receipt_hashes"]}
    assert len(expected) == ledger["parts"] and len(index["archives"]) == 4
    with tempfile.TemporaryDirectory(dir=a.scratch) as tmp:
        d = Path(tmp)
        for name in ("native", "deep"):
            (d / name).mkdir()
        for shard, obj in enumerate(index["archives"]):
            name = f"tomar-grouped-08371-{shard}.zip"
            path = a.bundle_index.parent / name
            assert obj["name"] == name and obj["sha256"] == sha(path.read_bytes())
            assert obj["bytes"] == path.stat().st_size < 10*1024**2
            entries = {}
            for part, hashes in expected.items():
                if part % 4 != shard:
                    continue
                stem = f"part{part:04}"
                for label, digest in zip([f"native/{stem}.json", f"native/{stem}.npz", f"deep/{stem}.json"], hashes):
                    entries[label] = digest
            with ZipFile(path) as archive:
                assert set(archive.namelist()) == set(entries) and len(archive.namelist()) == len(entries) == obj["members"]
                for label, digest in entries.items():
                    data = archive.read(label)
                    assert sha(data) == digest
                    (d / label).write_bytes(data)
        source = Path(__file__).with_name("audit_tomar_grouped_slice.py")
        replay = subprocess.check_output([sys.executable, "-I", "-B", source, a.certdir, a.table,
            a.slice_index, d/"native", d/"deep", d, a.cache, *(["--partial"] if a.partial else [])])
        assert json.loads(replay) == ledger
    print(json.dumps({"target": "0.8371", "complete": index["complete"], "parts": index["parts"],
          "bundle_index_sha256": sha(a.bundle_index.read_bytes()), "audit_sha256": index["audit_sha256"],
          "totals": ledger["totals"], "historical_area_inputs": [], "public_archives": index["archives"]},
          indent=2))
if __name__ == "__main__":
    main()
