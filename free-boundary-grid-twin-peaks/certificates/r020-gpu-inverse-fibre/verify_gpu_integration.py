#!/usr/bin/env python3
"""Exact-token regression and independent host-DP checks for the terminal GPU kernel."""
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent
OLD = ["tau", "replicate", "seed_hex", "Z", "U", "logZ", "U_over_Z", "scaled_U",
       "D", "F", "G", "F_over_Z", "G_over_Z", "V", "C", "F_minus_2G_minus_4D", "dead"]

def line(path):
    return next(json.loads(x) for x in path.read_text().splitlines() if '"kind":"replicate"' in x)

def token(text, name):
    m = re.search(r'"' + re.escape(name) + r'":(?:("(?:[^"\\]|\\.)*")|([^,}]+))', text)
    assert m, (name, text)
    return m.group(1) or m.group(2)

def exact_tokens(a, b):
    aa = next(x for x in a.read_text().splitlines() if '"kind":"replicate"' in x)
    bb = next(x for x in b.read_text().splitlines() if '"kind":"replicate"' in x)
    for k in OLD:
        assert token(aa, k) == token(bb, k), (a.name, b.name, k, token(aa, k), token(bb, k))

def host(path):
    # This is deliberately a separate Python translation of the stated DP.
    import subprocess
    subprocess.run(["python3", "host_inverse_from_trace.py", str(path) + ".trace.json"], cwd=ROOT,
                   check=True, stdout=subprocess.DEVNULL)
    return json.loads(Path(str(path) + ".inverse-host.json").read_text())

def check(stem):
    gpu = ROOT / (stem + ".jsonl")
    base = ROOT / ("base-" + stem + ".jsonl")
    exact_tokens(gpu, base)
    g, h = line(gpu), host(gpu)
    assert [str(x) for x in g["inverse_m2_groups"]] == h["m2_groups"], (stem, g, h)
    assert [str(x) for x in g["inverse_m3_groups"]] == h["m3_groups"], (stem, g, h)
    assert str(sum(g["inverse_m2_groups"])) == h["m2_total"]
    assert str(sum(g["inverse_m3_groups"])) == h["m3_total"]
    assert g["K2"] is None and h["overflow"] is False
    return {"case": stem, "groups": len(g["inverse_m2_groups"]), "m2_total": h["m2_total"], "m3_total": h["m3_total"],
            "mean_m2": g["mean_m2"], "mean_m3": g["mean_m3"]}

if __name__ == "__main__":
    cases = ["n3-p17-t10-l0", "n3-p17-t10-l4", "n3-p17-t10-l16", "n8-p257-t65-l4"]
    rows = [check(x) for x in cases]
    out = {"status": "passed", "base_numeric_tokens": "identical for all pre-existing replicate fields",
           "independent_python_dp": "per-terminal-group exact equality", "rows": rows,
           "uint64_bound": {"n_max": 16, "N_max": 256,
                             "per_word_m2": "4*(N-1) <= 1020",
                             "per_word_m3": "4*(N-1)*(N-2) <= 259080",
                             "terminal_weight": "a_i in {0,1}",
                             "per_group_m3": "256*259080 = 66324480", "less_than": "2^64"}}
    (ROOT / "verification.json").write_text(json.dumps(out, indent=2) + "\n")
    print(json.dumps(out, indent=2))
