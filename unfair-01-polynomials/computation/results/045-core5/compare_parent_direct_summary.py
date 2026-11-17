#!/Users/zzh/.venvs/main/bin/python
"""Compare already-produced summaries only; never imports either search code."""
from __future__ import annotations
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
OWN = HERE / "hgq-depth64-cap100000" / "summary.json"
PARENT = Path("/Volumes/data/CodexTemp/T2/unfair-01-polynomials/round045-parent-core5/direct-prefix-summary.json")

def sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()

def main() -> None:
    own, parent = json.loads(OWN.read_text()), json.loads(PARENT.read_text())
    assert len(own["results"]) == len(parent["rows"]) == 404
    for number, (lhs, rhs) in enumerate(zip(own["results"], parent["rows"])):
        assert lhs["core"] == rhs["core"], (number, lhs["core"], rhs["core"])
        assert lhs["status"] == "dead_exhaustive_prefix"
        assert lhs["death"]["depth"] == rhs["depth"], number
        assert [x["accepted"] for x in lhs["layers"]] == rhs["counts"], number
    result = {
        "status": "passed",
        "comparison": "all 404 cores: core tuple, every layer state count, and first death depth agree",
        "own_summary_sha256": sha(OWN),
        "parent_summary_sha256": sha(PARENT),
    }
    target = HERE / "hgq_vs_parent_comparison.json"
    target.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))

if __name__ == "__main__":
    main()
