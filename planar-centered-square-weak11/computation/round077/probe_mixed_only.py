#!/usr/bin/env python3
"""One finite FLOAT screen of the uncapped two-outside mixed-only candidate."""
import argparse
import concurrent.futures as cf
from fractions import Fraction as F
import hashlib
import importlib.util
import json
import math
import multiprocessing as mp
import os
from pathlib import Path
import platform
import random
import time


def init(path):
    global geom
    sp = importlib.util.spec_from_file_location('geom', path)
    geom = importlib.util.module_from_spec(sp)
    sp.loader.exec_module(geom)


def one(seed):
    rng = random.Random(seed)
    mass = [F(1)] + [F(rng.choice((1, 2, 4, 8, 16, 32, 64)),
                           rng.choice((1, 2, 4, 8, 16))) for _ in range(2)]
    reach = max(10, math.ceil(30 * math.sqrt(float(sum(mass)))))
    while True:
        pts = [(0, 0)]
        if seed % 4 == 0:
            pts += [(-rng.randint(1, reach), rng.randint(-5, 5)),
                    (rng.randint(-5, 5), -rng.randint(1, reach))]
        else:
            pts += [(rng.randint(-reach, reach), rng.randint(-reach, reach))
                    for _ in range(2)]
        if len(set(pts)) == 3:
            break
    p = [tuple(F(v, 20) for v in q) for q in pts]
    cap = math.sqrt(float(sum(mass))) / 2
    def receiver(ids):
        m = sum(mass[i] for i in ids)
        aa = [min(p[i][j] for i in ids) for j in (0, 1)]
        bb = [max(p[i][j] for i in ids) for j in (0, 1)]
        if max(bb[j] - aa[j] for j in (0, 1))**2 >= m:
            return None
        r = math.sqrt(float(m)) / 2
        x0, x1 = max(-cap, float(bb[0])-r), min(0, float(aa[0])+r)
        y0, y1 = max(-cap, float(bb[1])-r), min(0, float(aa[1])+r)
        return (x0, x1, y0, y1) if x0 < x1 and y0 < y1 else None
    full = [r for ids in ((0, 1), (0, 2), (0, 1, 2))
            if (r := receiver(ids)) is not None]
    residual = [r for ids in ((1,), (2,), (1, 2))
                if (r := receiver(ids)) is not None]
    area = geom.sweep_areas(full, residual)[2] if full else 0.0
    return {'seed': seed, 'positions': [list(map(str, q)) for q in p],
            'masses': list(map(str, mass)), 'mixed_area_float': area,
            'pid': os.getpid()}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--geometry', required=True)
    ap.add_argument('--output', type=Path, required=True)
    ap.add_argument('--count', type=int, default=1536)
    args = ap.parse_args()
    assert os.path.ismount('/Volumes/data')
    t = time.monotonic()
    with cf.ProcessPoolExecutor(3, mp_context=mp.get_context('spawn'),
                                initializer=init, initargs=(args.geometry,)) as pool:
        rows = list(pool.map(one, range(770000, 770000+args.count), chunksize=16))
    rows.sort(key=lambda r: r['mixed_area_float'], reverse=True)
    out = {'schema': 'r077-uncapped-mixed-only-float-screen',
           'python': platform.python_version(), 'count': len(rows),
           'parent_pid': os.getpid(), 'worker_pids': sorted({r['pid'] for r in rows}),
           'pool_joined': True, 'seed_start': 770000,
           'source_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
           'geometry_sha256': hashlib.sha256(Path(args.geometry).read_bytes()).hexdigest(),
           'top20': rows[:20],
           'count_above_half_plus_1e_minus_5': sum(r['mixed_area_float'] > .50001 for r in rows),
           'seconds': time.monotonic()-t,
           'scope': 'Rational input and exact eligibility; all areas FLOAT. Only full labels containing old and at least one outside, minus the complete residual. No old singleton. No grid/cap assumption. Finite screen is not a theorem.'}
    args.output.write_text(json.dumps(out, indent=2)+'\n')
    print(json.dumps({'count': len(rows), 'top': rows[0],
                      'count_above_half_plus_1e_minus_5': out['count_above_half_plus_1e_minus_5'],
                      'worker_pids': out['worker_pids'], 'seconds': out['seconds']}))


if __name__ == '__main__':
    main()
