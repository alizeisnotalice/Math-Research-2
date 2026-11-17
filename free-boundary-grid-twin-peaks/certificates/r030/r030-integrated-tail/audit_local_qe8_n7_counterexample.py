#!/usr/bin/env python3
"""Read the certified n=7 centre short-tail audit and check local QE8 failure."""
from __future__ import annotations
import json
from fractions import Fraction
from pathlib import Path


def main(source: Path, output: Path):
    x = json.loads(source.read_text())
    assert x["status"] == "passed" and x["n"] == 7 and x["v"] == 24
    T4 = int(x["residual_L_at_least_4_mass"])
    B = int(x["reference_selected_normalizer"]["Bv"])
    ratio, target = Fraction(T4, B), Fraction(8, 4 * 5)
    assert ratio > target
    output.write_text(json.dumps({"status": "passed", "scope": "exact finite n=7 centre only",
                                  "T_v_4_over_Bv": str(ratio), "target_8_over_4_times_5": str(target),
                                  "local_QE8_fails": True}, indent=2) + "\n")


if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser(); p.add_argument("source", type=Path); p.add_argument("output", type=Path)
    a = p.parse_args(); main(a.source, a.output)
