#!/usr/bin/env python3
"""Fresh all-profile n=5 gate for curved_single_region_root."""
import hashlib
import json
from pathlib import Path
import subprocess
import sys

binary = Path(sys.argv[1]).resolve()
reference = Path(sys.argv[2]).resolve()
out = Path(sys.argv[3]).resolve()
if out.exists():
    raise SystemExit("gate output directory must be fresh")
out.mkdir(parents=True)

def sha(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()

ref = json.loads(reference.read_text())
rows = {r["curve"]: r for r in ref["rows"]}
assert len(rows) == 32 and set(rows) == {"".join(x) for x in __import__("itertools").product("23", repeat=5)}
receipts = []
for curve in sorted(rows):
    target = out / f"{curve}.json"
    p = subprocess.Popen([str(binary), "5", str(target), curve])
    rc = p.wait()
    got = json.loads(target.read_text()) if rc == 0 and target.is_file() else None
    want = rows[curve]
    ok = (got is not None and got["curve"] == curve and
          got["left_vertices"] == want["left_vertices"] and
          got["left_rooted_count"] == want["left_rooted_count"] and
          got["right_rooted_count"] == want["right_rooted_count"] and
          got["count"] == want["count"])
    receipts.append({"curve": curve, "pid": p.pid, "returncode": rc, "waited_reaped": True,
                     "output_sha256": sha(target) if target.exists() else None, "matched": ok})
    if not ok:
        raise SystemExit(f"n5 mismatch for {curve}")
result = {"schema": "curved-single-region-n5-gate-v1", "profiles": len(receipts),
          "matched": all(r["matched"] for r in receipts), "binary_sha256": sha(binary),
          "reference_sha256": sha(reference), "receipts": receipts}
(out / "audit.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
print(json.dumps({k: v for k, v in result.items() if k != "receipts"}, sort_keys=True))
