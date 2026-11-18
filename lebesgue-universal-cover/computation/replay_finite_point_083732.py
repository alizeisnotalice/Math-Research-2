"""Replay the P30 terminal certificate on every R141 .83732 original leaf."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import ctypes
from fractions import Fraction as F
import json
from pathlib import Path
import platform
import shutil
import subprocess

import flint
from flint import arb, ctx
import numpy as np

from certify_finite_point_leaf import point_tables, upper_q
from certify_nine import ar
from certify_tomar_six_high import (
    FINE, Q, U, SourceTree, DirectionProof, build_kernel, ceil, digest, load_table,
)


HERE = Path(__file__).resolve().parent


class PointProof(DirectionProof):
    def __init__(self, cloud, cells, expected_diameter):
        ctx.prec = 160
        self.pi = arb.pi()
        assert len(cloud["points"]) == 30
        d2, ptr, dual = point_tables({**cloud, "cells": cells})
        assert d2 == expected_diameter and len(ptr) == 385 and ptr[-1] == 2493
        nums = cloud["directions"]
        self.v = [(ar(F(x, U)), ar(F(y, U)), ar(F(x*x+y*y, U*U)).sqrt()) for x,y in nums]
        self.directions = np.asarray(nums, dtype=np.int64).ravel()
        self.disk = np.asarray([upper_q(v[2]/2) for v in self.v], dtype=np.int64)
        self.pointers = np.asarray(ptr, dtype=np.int64)
        self.triples = np.asarray(dual, dtype=np.int64)


def bound_input(name, expected):
    raw = (HERE/name).read_bytes()
    assert digest(raw) == expected
    return json.loads(raw)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for k in ("tree","table32","table64","work","cache"):
        p.add_argument(k, type=Path)
    p.add_argument("--workers", type=int, default=3)
    p.add_argument("--parts", nargs="+", type=int)
    a = p.parse_args()

    data = (HERE/"finite-point-lift-083732.json").read_bytes()
    rec = json.loads(data)
    goal = F(rec["target"])
    assert goal == F(".83732")
    roots = bound_input("tomar-083732-retained.json", rec["retained_sha256"])
    audit = bound_input("tomar-083732-score-layers-audit.json", rec["localization_audit_sha256"])
    assert F(audit["target"]) == goal and audit["retained"] == roots
    old = bound_input("finite-point-lift-083731.json", rec["previous_614_sha256"])
    assert F(old["target"]) == F(".83731")
    known = {tuple(row[:2]) for row in old["rows"]}
    expected = sorted((int(part), pos) for part, ids in roots.items() for pos in ids)
    rows = rec["rows"]
    assert [tuple(r[:2]) for r in rows] == expected
    assert len(expected) == len(set(expected)) == 705 and len(known) == 614 and known <= set(expected)

    cloud = bound_input("finite-point-leaf-input.json", rec["point_input_sha256"])
    bank = bound_input("finite-point-global-input.json", rec["bank_sha256"])
    base = bound_input("tomar-six-high-input.json", rec["source_tables_sha256"])
    ctx.prec = 160
    pi = arb.pi()
    r = ar(F(707, 1000))
    assert r*r < ar(F(1, 2))
    assert (pi-(1/(2*r)).acos())/4+(r*r-ar(F(1, 4))).sqrt()/2 > ar(F(".838"))
    assert goal < F(".838")
    tree = SourceTree(a.tree)
    tbl = {n:load_table(getattr(a, f"table{n}"), tree, base["tables"][str(n)]) for n in (32,64)}
    point = PointProof(cloud, bank["cells"], rec["diameter_squared_numerator"])

    out = a.work.resolve()
    cache = a.cache.resolve()
    _, compiler_version = build_kernel(out, base, cache)
    src = HERE/"tomar_six_083730_guard.cpp"
    assert digest(src.read_bytes()) == rec["weighted_guard_sha256"]
    clang = shutil.which("clang++")
    assert clang
    library = out/"p30-lift-083732.so"
    subprocess.run([clang,"-O3","-std=c++17","-ftrapv","-shared","-fPIC",
        f"-fmodules-cache-path={cache}","-I",str(out),str(src),"-o",str(library)],check=True)
    guard = ctypes.CDLL(str(library)).lift_guard_run
    vector = lambda d:np.ctypeslib.ndpointer(dtype=d, ndim=1, flags="C_CONTIGUOUS")
    i64, i32 = vector(np.int64), vector(np.int32)
    guard.argtypes = [i64,i64,i64,i32,i32,ctypes.c_int64,ctypes.c_int,ctypes.c_int,
                      i64,i64,i64,i64,i64,ctypes.c_int,ctypes.c_int,ctypes.c_int,i64]
    guard.restype = ctypes.c_int

    def scaled(values):
        x = [F(v)*FINE for v in values]
        assert all(q.denominator == 1 for q in x)
        return np.array(list(map(int,x)), dtype=np.int64)

    root = scaled(tree.root["hi"])
    target = goal*2*Q*Q
    assert target.denominator > 1 and ceil(target) == 120670529279675639

    def run(row, box, support):
        part,pos,n,calls,old_count,new_count,depth,attempted,joint_count = row
        lo,hi = (scaled(x) for x in zip(*box,strict=True))
        status = np.zeros(30,dtype=np.int64)
        assert guard(lo,hi,root,*tbl[n],ceil(target),depth+1,calls,*support,
                     len(support[1]),len(support[3])-1,2,status) == 0
        assert status[:7].tolist() == [1,calls,old_count,new_count,depth,0,attempted], (part,pos,status[:7])
        assert calls == 2*(old_count+new_count)-1 and int(status[29]) == joint_count
        return row

    active = [x for x in rows if not a.parts or x[0] in a.parts]
    assert active
    results = []
    with ThreadPoolExecutor(a.workers) as pool:
        pending = []
        for row in active:
            box = tree.leaf(*row[:2])
            pending.append(pool.submit(run,row,box,point.presolve(box)))
            if len(pending) >= a.workers*2:
                results.append(pending.pop(0).result())
        results.extend(t.result() for t in pending)
    assert results == active
    print(json.dumps({
        "threshold":rec["target"], "complete":len(results)==705,
        "covered":len(results), "from_prior_614":sum(tuple(r[:2]) in known for r in results),
        "new_exterior":sum(tuple(r[:2]) not in known for r in results),
        "P30_used_leaves":sum(r[5]>0 for r in results),
        "nodes":sum(r[3] for r in results), "old_terminals":sum(r[4] for r in results),
        "P30_terminals":sum(r[5] for r in results), "joint_cap_terminals":sum(r[8] for r in results),
        "input_sha256":digest(data), "retained_sha256":rec["retained_sha256"],
        "compiler":compiler_version,"python":platform.python_version(),"python_flint":flint.__version__,
    },indent=2))


if __name__ == "__main__":
    main()
