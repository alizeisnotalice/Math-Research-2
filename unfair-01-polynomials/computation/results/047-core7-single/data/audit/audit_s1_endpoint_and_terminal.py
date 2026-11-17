#!/usr/bin/env python3
"""Small exact audit for the degree-seven S={1} terminal HGQ certificate.

It independently counts the trace/norm input box and derives the final fifth
Newton sum in Z[a]/(a^2+2a-2).  It does not repeat the archived full prefix
enumeration; that enumeration is bound by its source and summary hashes.
"""

from __future__ import annotations

import hashlib
import json
from math import isqrt
from pathlib import Path


ROOT = Path("/Volumes/data/CodexTemp/T2/unfair-01-polynomials/round047-parent")


def add(x, y):
    return (x[0] + y[0], x[1] + y[1])


def scale(k, x):
    return (k * x[0], k * x[1])


def mul(x, y):
    # a^2=-2a+2, corresponding to T=N=-2 in a^2=T*a-N.
    u, v = x
    s, t = y
    return (u * s + 2 * v * t, u * t + v * s - 2 * v * t)


def power(a, n):
    out = (1, 0)
    for _ in range(n):
        out = mul(out, a)
    return out


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def inputs():
    for t in range(-16, 17):
        for n in range(-72, 73):
            delta = t * t - 4 * n
            if 0 < delta < (17 - abs(t)) ** 2 and isqrt(delta) ** 2 != delta:
                yield (t, n, delta)


def main():
    rows = list(inputs())
    assert len(rows) == 680
    assert (-2, -2, 12) in rows
    a = (0, 1)
    # s5=-a^5+5a^3*c2-5a*c2^2-5a^2*c3+5*c2*c3+5a*c4-5*c5.
    s5 = add(
        add(scale(-1, power(a, 5)), scale(5, power(a, 3))),
        add(scale(-5, a), scale(-20, power(a, 2))),
    )
    s5 = add(s5, (20 - 25, 5))
    assert s5 == (-33, 26)
    # Conjugate a=-1-sqrt(3), so s5bar=-59-26sqrt(3).
    # 59+26sqrt(3) > 3107/32 iff 832sqrt(3)>1219.
    assert 3 * 832**2 > 1219**2
    summary = json.loads((ROOT / "s1-hgq-all680-summary.json").read_text())
    last = json.loads((ROOT / "s1-last-core-p5.json").read_text())
    ledger = [json.loads(line) for line in (ROOT / "s1-hgq-all680.jsonl").read_text().splitlines()]
    assert len(ledger) == 680
    assert {tuple(row["TN"]) for row in ledger} == set(rows)
    aggregate = {key: sum(row["counts"][key] for row in ledger)
                 for key in ("c2", "c3", "c4", "prefix4", "c5", "prefix6")}
    assert summary["TN"] == 680
    assert aggregate == summary["counts"] == {"c2": 2191, "c3": 6969, "c4": 62117,
                                                "prefix4": 23673, "c5": 1, "prefix6": 1}
    terminal = [row for row in ledger if row["counts"]["c5"]]
    assert len(terminal) == 1 and terminal[0]["TN"] == [-2, -2, 12]
    assert last["power_sums"][-1] == [-33, 26]
    assert last["p5_bound"] == [3107, 32] and not last["passes_both_p5"]
    out = {
        "status": "passed_terminal_arithmetic_and_input_coverage",
        "TN_count": len(rows),
        "terminal_TND": [-2, -2, 12],
        "s5_in_basis_1_a": list(s5),
        "integer_witness": {"left": 3 * 832**2, "right": 1219**2},
        "parent_summary_sha256": sha(ROOT / "s1-hgq-all680-summary.json"),
        "parent_680_ledger_sha256": sha(ROOT / "s1-hgq-all680.jsonl"),
        "parent_terminal_sha256": sha(ROOT / "s1-last-core-p5.json"),
        "parent_verifier_sha256": sha(ROOT / "verify_s1_hgq.py"),
        "scope": "Does not rerun the 680-core HGQ prefix enumeration; checks its input box and terminal certificate.",
    }
    target = Path(__file__).with_name("audit_s1_endpoint_and_terminal_output.json")
    target.write_text(json.dumps(out, indent=2) + "\n")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
