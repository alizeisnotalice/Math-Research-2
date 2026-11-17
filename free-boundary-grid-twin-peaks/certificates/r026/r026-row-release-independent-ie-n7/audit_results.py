#!/usr/bin/env python3
"""Verify the downloaded finite selected-boundary IE receipts."""
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
R = HERE / "results"
EXPECTED_SOURCE = "850dea773006f6fccecaf3ded0980604071da056c01dd595ded5873a3f558f42"
EXPECTED_BRIDGE = {
    "n6-toprow-bridge.json": "20779014929894865426585600",
    "n6-tworows-bridge.json": "6637177983794054400",
}

def sha(p):
    h = hashlib.sha256()
    with p.open("rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()

def resource(text):
    out = {}
    for line in text.splitlines():
        if " real" in line:
            out["time_line"] = line.strip()
        if "maximum resident set size" in line:
            out["maximum_resident_set_size_bytes"] = int(line.split()[0])
        if line.strip().endswith("swaps"):
            out["swaps"] = int(line.split()[0])
        if "peak memory footprint" in line:
            out["peak_memory_footprint_bytes"] = int(line.split()[0])
    return out

audit = {
    "scope": "Strict finite exact selected-boundary inclusion-exclusion computations only; not a proof about asymptotic twin peaks.",
    "source_sha256_expected": EXPECTED_SOURCE,
    "source_sha256_downloaded": sha(HERE / "completion_ie_selected.cpp"),
    "n6_bridge": {},
    "n7": {},
}
for file, expected in EXPECTED_BRIDGE.items():
    d = json.loads((R / file).read_text())
    audit["n6_bridge"][file] = {"hS": d["hS"], "expected_hS": expected,
                                "matches": d["hS"] == expected, "candidate_holds": d["candidate_holds"]}
for tag in ("n7-toprow", "n7-tworows"):
    d = json.loads((R / (tag + ".json")).read_text())
    e = json.loads((R / (tag + ".exit.json")).read_text())
    rs = R / (tag + ".resource.txt")
    audit["n7"][tag] = {
        "result": d,
        "exit": e,
        "resource": resource(rs.read_text()),
        "downloaded_output_sha256": sha(R / (tag + ".json")),
        "output_hash_matches_receipt": sha(R / (tag + ".json")) == e["output_sha256"],
        "resource_hash_matches_receipt": sha(rs) == e["resource_sha256"],
    }
audit["all_checks_pass"] = (
    audit["source_sha256_downloaded"] == EXPECTED_SOURCE
    and all(v["matches"] for v in audit["n6_bridge"].values())
    and all(v["exit"]["exit_code"] == 0 and v["output_hash_matches_receipt"] and v["resource_hash_matches_receipt"]
            and v["resource"].get("swaps") == 0 for v in audit["n7"].values())
)
(HERE / "audit-receipt.json").write_text(json.dumps(audit, indent=2, sort_keys=True) + "\n")
print(json.dumps({"all_checks_pass": audit["all_checks_pass"], "n7": audit["n7"]}, indent=2))
