#!/usr/bin/env python3
import csv
import hashlib
import json
import math
import pathlib
import statistics

root = pathlib.Path(__file__).resolve().parent
old = list(csv.DictReader((root / "reference" / "original-reps0-47.csv").open()))
ledger = list(csv.DictReader((root / "canonical-64-ledger.csv").open()))
receipt = json.loads((root / "recovery-receipt.json").read_text())
base, vertex, L = 0xA024D0000000, 119, 16

def sha(p):
    return hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest()
def weights(rows):
    z = [float(x["probZ"]) for x in rows]
    s = sum(z); w = [x/s for x in z]
    se = statistics.stdev(z) / math.sqrt(len(z))
    return {"replicate_weight_ess": 1.0/sum(x*x for x in w), "max_share": max(w),
            "mean_probZ": s/len(z), "relative_sample_se": se/(s/len(z))}

assert len(old) == 48 and [int(x["rep"]) for x in sorted(old, key=lambda x:int(x["rep"]))] == list(range(48))
assert len(ledger) == 64 and [int(x["rep"]) for x in ledger] == list(range(64))
assert all(int(x["seed"]) == base + vertex*1000000007 + int(x["rep"])*1009 + L*1000003 for x in ledger)
assert all(int(x["dead"]) == 0 and math.isfinite(float(x["probZ"])) and float(x["probZ"]) > 0 for x in ledger)
assert all(x["terminal_hash"].isdigit() and int(x["terminal_hash"]) != 0 for x in ledger)
old_by_rep = {int(x["rep"]): x for x in old}
assert all(ledger[q][k] == old_by_rep[q][k] for q in range(48) for k in old_by_rep[q].keys())
for q in range(48, 64):
    raw = root / "raw" / f"recovered-target-rep{q:02d}.csv"
    row = list(csv.DictReader(raw.open()))
    assert len(row) == 1 and int(row[0]["rep"]) == 0
    assert ledger[q]["source_raw"] == raw.name and ledger[q]["source_local_rep"] == "0"
    assert all(ledger[q][k] == row[0][k] for k in row[0].keys() if k != "rep")
    assert int(row[0]["seed"]) == base + vertex*1000000007 + q*1009 + L*1000003
    resource = (root / "raw" / f"recovered-target-rep{q:02d}.resource.txt").read_text()
    assert "real" in resource and "maximum resident set size" in resource
assert receipt["returncode"] == 0 and not (root / "exception-cleanup.json").exists()
old_p65 = list(csv.DictReader((root / "reference" / "comparison-p65536.csv").open()))
report = {
    "strict": {"old_raw_48_preserved_in_ledger": True, "canonical_64_reps": True,
        "seed_formula_64": True, "dead_zero_64": True, "positive_finite_normalizer_64": True,
        "terminal_hash_present_nonzero_64": True, "recovered_raw_and_resource_16": True,
        "driver_returncode_zero": True, "no_exception_cleanup": True,
        "canonical_sha256_matches_receipt": sha(root / "canonical-64-ledger.csv") == receipt["canonical_sha256"]},
    "terminal_hash": {"unique_count": len({x["terminal_hash"] for x in ledger}),
                      "note": "hashes are execution diagnostics; source asserts one peak per terminal particle"},
    "P262144": weights(ledger), "P65536": weights(old_p65),
    "hashes": {"canonical": sha(root / "canonical-64-ledger.csv"),
               "original_interrupted_raw": sha(root / "reference" / "original-reps0-47.csv")},
}
print(json.dumps(report, indent=2, sort_keys=True))
