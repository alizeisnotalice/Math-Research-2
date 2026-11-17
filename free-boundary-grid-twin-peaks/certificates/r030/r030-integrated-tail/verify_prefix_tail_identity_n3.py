#!/usr/bin/env python3
"""Independent-set verification of the R030 prefix-tail identity on all n=3 holes."""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path


SHORT_IE = None


def load_short_ie():
    spec = importlib.util.spec_from_file_location("short_ie", SHORT_IE)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def tail_formula(n: int, v: int, max_t: int):
    s = load_short_ie()
    N, M, adj = n * n, n * n - 1, s.graph(n)
    H, Nv, Z, bags, C, cls = s.c_h_values(n, v, adj)
    ans = []
    for t in range(1, max_t + 1):
        value = 0
        den = s.falling(M, t)
        for I, count in C.items():
            U, P = H & ~cls[I], Z & (H & ~cls[I])
            a = sum((P & ~s.closure(eta, adj)).bit_count()
                    for eta in s.eta_words(t - 1, U, Nv, adj))
            if a:
                assert count % den == 0
                value += (-1 if I.bit_count() & 1 else 1) * (count // den) * a
        ans.append(value)
    return ans


def main(tail_path: Path, output: Path):
    x = json.loads(tail_path.read_text())
    assert x["n"] == 3
    rows = []
    for q in x["holes"]:
        hist = list(map(int, q["histogram_L"]))
        expected = [sum(hist[t - 1:]) for t in range(1, len(hist) + 1)]
        got = tail_formula(3, q["v"], len(hist))
        assert got == expected, (q["v"], got, expected)
        rows.append({"v": q["v"], "all_tail_levels_match": True})
    output.write_text(json.dumps({"status": "passed", "n": 3, "holes": rows}, indent=2) + "\n")


if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser(); p.add_argument("tail", type=Path); p.add_argument("output", type=Path); p.add_argument("--short-ie",type=Path,required=True)
    a = p.parse_args(); SHORT_IE=a.short_ie; main(a.tail, a.output)
