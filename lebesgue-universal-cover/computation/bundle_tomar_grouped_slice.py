"""Seal audited native scores and deep original-leaf checks as four portable ZIP parts."""
import argparse
import hashlib
import json
from pathlib import Path
from zipfile import ZIP_DEFLATED, ZipFile, ZipInfo

def sha(data):
    return hashlib.sha256(data).hexdigest()

def main():
    p = argparse.ArgumentParser(description=__doc__)
    for key in ("audit", "native", "deep", "out"):
        p.add_argument(key, type=Path)
    a = p.parse_args()
    report = json.loads(a.audit.read_bytes())
    assert report["historical_area_inputs"] == [] and report["target"] == "0.8371"
    a.out.mkdir(parents=True, exist_ok=True)
    parts = report["receipt_hashes"]
    assert len(set(x[0] for x in parts)) == len(parts) == report["parts"]
    objects = []
    for shard in range(4):
        final = a.out / f"tomar-grouped-08371-{shard}.zip"
        tmp = final.with_suffix(".partial")
        count = 0
        with ZipFile(tmp, "w") as dest:
            for part, h1, h2, h3 in parts:
                if part % 4 != shard:
                    continue
                key = f"part{part:04}"
                for label, path, expected in (("native", a.native / (key + ".json"), h1),
                        ("native", a.native / (key + ".npz"), h2), ("deep", a.deep / (key + ".json"), h3)):
                    data = path.read_bytes()
                    assert sha(data) == expected
                    if path.suffix == ".json":
                        assert all(x not in data for x in (b"/Users/", b"/Volumes/", b"/home/"))
                    entry = ZipInfo(label + "/" + path.name)
                    entry.external_attr = 0o100644 << 16
                    dest.writestr(entry, data, compress_type=ZIP_DEFLATED, compresslevel=9)
                    count += 1
        with ZipFile(tmp) as check:
            assert check.testzip() is None and len(check.namelist()) == count
        assert tmp.stat().st_size < 10 * 1024**2
        tmp.replace(final)
        objects.append({"name": final.name, "members": count, "bytes": final.stat().st_size,
                        "sha256": sha(final.read_bytes())})
    result = {"target": report["target"], "complete": report["complete"], "parts": report["parts"],
              "audit_sha256": sha(a.audit.read_bytes()), "archives": objects}
    output = a.out / "tomar-grouped-08371-index.json"
    output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result))
if __name__ == "__main__":
    main()
