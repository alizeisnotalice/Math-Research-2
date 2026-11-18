"""Rebuild the exact two-polygon cover result in a new directory."""
from argparse import ArgumentParser
from hashlib import sha256
import json
from pathlib import Path
import shutil
import subprocess
import sys
from zipfile import ZipFile

STEPS = (
    "certify", "fine_certify", "derivative", "coarse_derivative", "micro_certify",
    "micro_derivative", "triple", "local_exact", "algebraic", "join_exact",
)
CANDIDATES = {"scan.json", "fine-scan.json", "micro-scan.json", "selected.json", "micro-selected.json"}


def replay(work):
    project = Path(__file__).resolve().parent
    inputs = json.loads((project / "inputs.json").read_bytes())
    for label in ("P", "Q"):
        path = project.parent / inputs[f"source_{label}"]
        source = path.read_bytes()
        assert sha256(source).hexdigest() == inputs[f"source_{label}_sha256"]
        original = json.loads(source)
        assert inputs[label] == (original[-1] if label == "P" else original)
    work.mkdir(parents=True, exist_ok=False)
    shutil.copyfile(project / "inputs.json", work / "inputs.json")
    hashes = json.loads((project / "candidate-hashes.json").read_bytes())
    with ZipFile(project / "candidates.zip") as archive:
        assert set(hashes) == set(archive.namelist()) == CANDIDATES
        assert len(archive.namelist()) == len(CANDIDATES)
        for name in sorted(CANDIDATES):
            data = archive.read(name)
            assert sha256(data).hexdigest() == hashes[name]
            (work / name).write_bytes(data)
    for step in STEPS:
        script = work / (step + ".py")
        shutil.copyfile(project / script.name, script)
        with (work / (step + ".log")).open("wb") as log:
            subprocess.run([sys.executable, "-I", "-B", str(script)], check=True, stdout=log)
    for name in ("triple.json", "local-exact.json", "joined-exact.json", "algebraic.json"):
        assert json.loads((work / name).read_bytes()) == json.loads((project / name).read_bytes())
    subprocess.run([sys.executable, "-I", "-B", str(project / "exact_hull.py"), str(work)], check=True)
    print(json.loads((work / "algebraic.json").read_bytes())["area_interval"])


if __name__ == "__main__":
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("work", type=Path)
    replay(parser.parse_args().work.resolve())
