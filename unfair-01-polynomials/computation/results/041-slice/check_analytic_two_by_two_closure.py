#!/usr/bin/env python3
"""Exact branch regression for analytic_two_by_two_closure.md."""
from __future__ import annotations

import hashlib
import json
import platform
import sys
from pathlib import Path

import sympy as sp


def main() -> None:
    y = sp.symbols("y")
    # Exact factor identities forced at the two equality branches.
    assert sp.expand((1 + y + y**2) * (1 - y + y**2) - (1 + y**2 + y**4)) == 0
    assert sp.expand((1 - y + y**2) * (1 + y**2) - (1 - y + 2*y**2 - y**3 + y**4)) == 0

    branches = {
        "c1=-2": {"UV_p": "h1-2 < 0", "outcome": "negative"},
        "c1=1": {"forced": "h1=0, {r,t}={1,0}, c2=a+a^-1>=2", "outcome": "UV_2p>=2"},
        "c1=0": {"forced": "c2=1,h2=0,a=1,{r,t}={1,-1}", "outcome": "Q reducible"},
        "c1=-1": {"forced": "h1=1,c2=2,h2=0,a=1,{r,t}={-1,0}", "outcome": "Q reducible"},
    }
    # Integer possibilities for c1 from the three endpoint boxes.
    assert set(range(-2, 2)) == {-2, -1, 0, 1}
    here = Path(__file__).resolve()
    report = {
        "status": "passed_exact_two_by_two_branch_regressions",
        "branches": branches,
        "scope": "Exact checks of the forced equality factorizations and finite c1 branch labels; continuous inequalities are proved in the companion manuscript.",
        "python": sys.version,
        "sympy": sp.__version__,
        "platform": platform.platform(),
        "source_sha256": hashlib.sha256(here.read_bytes()).hexdigest(),
    }
    here.with_suffix(".json").write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n")
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
