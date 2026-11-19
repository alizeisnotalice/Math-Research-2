"""Bind the cross-test to the physical upper construction and the fixed five-curve pose."""
import argparse
from fractions import Fraction
import json
from pathlib import Path
import subprocess
import sys
from zipfile import ZipFile

from flint import arb

ROOT = Path(__file__).resolve().parent
COMP = ROOT.parent


def load_json(path):
    return json.loads(path.read_bytes())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("work", type=Path)
    args = parser.parse_args()
    work = args.work.resolve()
    work.mkdir(parents=True, exist_ok=False)
    with ZipFile(ROOT / "new-spur.zip") as archive:
        archive.extractall(work / "spur-source")
    subprocess.run([sys.executable,"-I","-B",work/"spur-source/kit/replay.py",COMP/"released-gibbs-sector/base.zip",work/"physical"],check=True)
    geom = work / "physical/spur-proof/witness.py"
    variables = {"__file__":str(geom),"__name__":"frozen_five_vertices"}
    prefix, marker, _ = geom.read_text().partition("w,e=")
    assert marker
    exec(compile(prefix,str(geom),"exec"),variables)
    A = variables["A"]
    vertices,_ = variables["body"](variables["env"], A(Fraction(241,400)), A(-Fraction(22,1125)))
    with ZipFile(ROOT / "cross-pentagon.zip") as archive:
        archive.extractall(work/"cross")
    with ZipFile(COMP/"directional-seven/six-low-continuous.zip") as archive:
        match = [x for x in archive.namelist() if x.split("/")[-1]=="true-best-pose.json"]
        assert len(match)==1
        assert archive.read(match[0]) == (work/"cross/true-best-pose.json").read_bytes()
    assert (work/"cross/finite-point-leaf-input.json").read_bytes()==(COMP/"finite-point-leaf-input.json").read_bytes()
    frozen = load_json(work/"cross/spur-vertices.json")
    cloud = load_json(ROOT/"five-inner-points.json")
    assert len(vertices)==len(frozen)==len(cloud["points"])==5
    for pt,box in zip(vertices,frozen):
        assert all(arb(s).contains(x) for x,s in zip(pt,box))
    den = cloud["denominator"]
    for x,y in cloud["points"]:
        for a,b in vertices:
            assert 1-(arb(x)/den-a)**2-(arb(y)/den-b)**2>arb("9.4389e-9")
    record = (work/"cross/exclusion.json").read_bytes()
    subprocess.run([sys.executable,"-I","-B","certify_exclusion.py"],cwd=work/"cross",check=True)
    assert (work/"cross/exclusion.json").read_bytes()==record
    with ZipFile(ROOT/"illumination.zip") as archive:
        archive.extractall(work/"illumination")
    for source in sorted((work/"illumination").glob("plan-*.json")):
        stem=source.stem.removeprefix("plan-")
        out=work/("illumination-"+stem+".json")
        subprocess.run([sys.executable,"-B",ROOT/"verify_illumination.py",COMP,work/"cross/true-best-pose.json",ROOT/"five-inner-points.json",source,out],check=True)
        assert load_json(out)==load_json(work/"illumination"/("receipt-"+stem+".json"))


if __name__ == "__main__":
    main()
