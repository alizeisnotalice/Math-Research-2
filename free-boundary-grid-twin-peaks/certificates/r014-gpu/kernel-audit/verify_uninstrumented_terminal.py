#!/usr/bin/env python3
"""Check uninstrumented production CLI terminals against the rank audit."""
import json
from pathlib import Path

root = Path(__file__).parent
audit = json.loads((root / "audit-results.json").read_text())
rows = []
for n in (9, 10, 12):
    needle = f"n{n} tau2 P17 strat"
    expected = next(x for x in audit["records"] if x["case"] == needle)
    raw = [json.loads(line) for line in (root / f"original-n{n}-tau2-p17-strat-domain1.jsonl").read_text().splitlines()]
    rep = next(x for x in raw if x["kind"] == "replicate")
    assert raw[-1]["kind"] == "run_end"
    assert rep["seed_hex"] == expected["seed_hex"]
    assert rep["dead"] == expected["dead"] is False
    assert rep["Z"] == expected["Z"]
    assert rep["U"] == expected["U"]
    rows.append({"n": n, "seed_hex": rep["seed_hex"], "Z": rep["Z"], "U": rep["U"]})

result = {"status": "passed", "unmodified_cli_cases": 3, "rows": rows}
(root / "original-terminal-check.json").write_text(json.dumps(result, indent=2) + "\n")
print(json.dumps(result, indent=2))
