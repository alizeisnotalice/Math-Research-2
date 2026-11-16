#!/usr/bin/env python3
"""Export every n=3 one-peak word's exact m2,m3 from the reference DP."""
import runpy
from pathlib import Path


REFERENCE = str(Path(__file__).resolve().parents[1] / "r020-inverse-fibre" / "exact_inverse_fibre_n3.py")
g = runpy.run_path(REFERENCE)
rows = []
for word in g["onepeaks"]:
    rows.append(("-".join(map(str, word)),
                 g["inverse_dp"](word, 2) if g["geo_prefix"](word, 2) else 0,
                 g["inverse_dp"](word, 3) if g["geo_prefix"](word, 3) else 0))
assert len(rows) == 29568
path = Path(__file__).with_name("n3-python-fibres.csv")
path.write_text("word,m2,m3\n" + "".join(f"{w},{m2},{m3}\n" for w, m2, m3 in rows))
print(path, len(rows))
