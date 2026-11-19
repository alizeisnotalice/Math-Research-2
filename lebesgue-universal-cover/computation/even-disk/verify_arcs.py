"""Bind the two finite-arc completions to their normalized R275 point groups."""
import argparse
from fractions import Fraction
import json
from pathlib import Path
import subprocess
import sys
from zipfile import ZipFile

here = Path(__file__).resolve().parent
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("work", type=Path)
work = parser.parse_args().work.resolve()
work.mkdir(parents=True, exist_ok=False)
with ZipFile(here / "arc-completions.zip") as packed:
    packed.extractall(work)
with ZipFile(here.parent / "directional-seven" / "scale-200.zip") as source:
    original = source.read("points.json")
if original != (work / "points-source.json").read_bytes():
    raise ValueError("the archived completion has different source points")
sources = json.loads(original)
scale = Fraction(4999999, 5000000)
for key, count in (("P66",25), ("physicalQ",11)):
    group = json.loads((work / (key+".json")).read_bytes())
    fixed = sources[key]
    assert group["points"] == fixed["points"]
    assert Fraction(1,group["denominator"]) == Fraction(1,fixed["denominator"])/scale
    subprocess.run([sys.executable,"-I","-B","verify_completion.py",key],cwd=work,check=True)
    result = json.loads((work / (key+"-completion.json")).read_bytes())
    assert result["complete"] and result["completions"][0]["unit_width_arcs"] == count
report = subprocess.check_output([sys.executable,"-I","-B","thresholds.py"],cwd=work)
result = json.loads(report)["rows"]
frozen = json.loads((work / "thresholds.json").read_bytes())["rows"]
assert [(x["method"],x["target"],x["first_odd"]) for x in result] == [(x["method"],x["target"],x["first_odd"]) for x in frozen]
(work / "new-thresholds.json").write_bytes(report)
print("Verified both unit-width completions and all threshold inequalities")
