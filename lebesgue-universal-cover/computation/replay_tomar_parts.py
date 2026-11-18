"""Replay selected Tomar v9 parts and their exact trees from a verified archive.

Usage: python -B replay_tomar_parts.py UPSTREAM EXTRACTED_CERTDIR [0,1,...]
Omit the part list to verify the full certificate. stdout contains a JSON receipt.
python-flint and numpy must be installed in this interpreter.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
import tarfile
import time

import flint
import numpy as np


def digest(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("upstream", type=Path)
    parser.add_argument("certdir", type=Path)
    parser.add_argument("parts", nargs="?")
    args = parser.parse_args()
    source = args.upstream.resolve()
    directory = args.certdir.resolve()
    pinned = "10af98606c98ca9ffd0d9f402cae1e2efe74b792"
    commit = subprocess.check_output(["git", "-C", str(source), "rev-parse", "HEAD"], text=True).strip()
    if commit != pinned:
        raise ValueError("unexpected upstream commit")
    archive = directory.with_name(directory.name + ".tar.xz")
    manifests = source / "certificates"
    if not archive.is_file():
        archive = manifests / archive.name
    expected = []
    for manifest in [manifests/"SHA256SUMS", manifests/"cert11d_08366.sha256"]:
        for line in manifest.read_text().splitlines():
            value, name = line.split(maxsplit=1)
            if name.strip() == archive.name:
                expected.append(value)
    if len(expected) != 1 or digest(archive) != expected[0]:
        raise ValueError("unexpected archive checksum")
    with tarfile.open(archive, "r|xz") as contents:
        for item in contents:
            if item.name.removeprefix("./") == directory.name + "/SHA256SUMS":
                embedded = contents.extractfile(item).read()
                break
        else:
            raise ValueError("archive lacks manifest")
    if (directory/"SHA256SUMS").read_bytes() != embedded:
        raise ValueError("extracted manifest differs from pinned archive")
    files = embedded.decode().splitlines()
    for line in files:
        value, name = line.split(maxsplit=1)
        if digest(directory/name.strip()) != value:
            raise ValueError("incorrect checksum for " + name)
    scripts = {name: source/"checker"/name for name in ("treecheck.py", "lbcheck.py")}
    for name, path in scripts.items():
        raw = subprocess.check_output(["git","-C",str(source),"show",f"{pinned}:checker/{name}"])
        if hashlib.sha256(raw).hexdigest() != digest(path):
            raise ValueError("modified upstream checker " + name)

    headers = sorted(directory.glob("part*.json"))
    header = json.loads(headers[0].read_text())
    n = header["nparts"]
    parts = list(range(n)) if args.parts is None else sorted(set(map(int,args.parts.split(","))))
    if not parts or parts[0] < 0 or parts[-1] >= n:
        raise ValueError("part out of range")
    env = {**os.environ, "OPENBLAS_NUM_THREADS":"1", "OMP_NUM_THREADS":"1", "VECLIB_MAXIMUM_THREADS":"1"}

    def run(name, *extra):
        command = [sys.executable, "-I", "-B", str(scripts[name]), str(directory), *extra]
        t = time.monotonic()
        child = subprocess.run(command, text=True, capture_output=True, check=True, env=env)
        return round(time.monotonic()-t,3), [json.loads(x) for x in child.stdout.splitlines()]

    initial_time, presence = run("lbcheck.py","present")
    if not (len(presence)==1 and presence[0][0] is True):
        raise ValueError("incomplete metadata")
    trees, balls = [], []
    for p in parts:
        t, lines = run("treecheck.py",str(p),str(p+1))
        tree = lines[0]
        if (len(lines)!=1 or tree["ok"] is not True or tree["first"]!=p or tree["last"]!=p+1):
            raise ValueError("tree failure at " + str(p))
        t_ball, lines = run("lbcheck.py",str(p))
        ball = lines[0]
        if (len(lines)!=1 or ball["ok"] is not True or ball["fails"]!=0 or
                ball["part"]!=p or ball["leaves"]!=tree["leaves"]):
            raise ValueError("ball/partition mismatch at " + str(p))
        trees.append({"seconds":t,**{k:tree[k] for k in ("cells","leaves","fan_leaves","far_leaves","multi_leaves","max_depth")}})
        balls.append({"seconds":t_ball,**ball})
    full = parts == list(range(n))
    print(json.dumps({"upstream":commit,"python":platform.python_version(),"numpy":np.__version__,
        "python_flint":flint.__version__,"archive_sha256":expected[0],"verified_file_count":len(files),
        "checker_sha256":{name:digest(path) for name,path in scripts.items()},
        "header":header,"presence":presence[0],"presence_seconds":initial_time,
        "complete":full,"parts":parts,"total_leaves":sum(t["leaves"] for t in trees),
        "trees":trees,"balls":balls},indent=2))

if __name__ == "__main__":
    main()
