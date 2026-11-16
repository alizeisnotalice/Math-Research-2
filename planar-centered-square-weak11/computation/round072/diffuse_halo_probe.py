#!/usr/bin/env python3
r"""Deterministic floating midpoint probe, NOT an area certificate.

S_n is the equal n^2-atom grid; nu=t*1_(B\J)dx with t=n/(n+1),
J=[-1/n^2,1]^2 and B=[-1,2]^2. For relevant centres and r<1/2,
Q is contained in B. The exact square-cap and finite-event arguments
are recorded separately. Only point samples and floating arithmetic run here.
"""
from __future__ import annotations
import argparse
import concurrent.futures
import hashlib
import json
import os
from pathlib import Path
import platform
import sys
import time
import numpy as np


def probe(task):
    n, resolution = task
    started = time.monotonic()
    w = 1.0-1.0/n**2
    h = (n+1.0)/n**2
    t = n/(n+1.0)
    jlo = -1.0/n**2
    grid = h*np.arange(n)
    axis = -.5+(np.arange(resolution)+.5)*(w+1)/resolution
    xx, yy = np.meshgrid(axis, axis, indexing='ij')
    points = np.column_stack([xx.ravel(), yy.ravel()])
    hit = internal = near = 0
    global_margin = -np.inf
    for start in range(0, len(points), 512):
        q = points[start:start+512]
        x, y = q[:, 0, None], q[:, 1, None]
        radii = np.concatenate((np.abs(x-grid), np.abs(y-grid)), axis=1)
        counts = []
        for center in (x, y):
            lower = np.maximum(0, np.ceil((center-radii)/h-1e-11))
            upper = np.minimum(n-1, np.floor((center+radii)/h+1e-11))
            counts.append(np.maximum(0, upper-lower+1))
        mass = counts[0]*counts[1]/n**2
        overlaps = []
        for center in (x, y):
            overlaps.append(np.maximum(0, np.minimum(center+radii, 1)
                                       -np.maximum(center-radii, jlo)))
        cost = t*overlaps[0]*overlaps[1]+(1-t)*4*radii**2
        valid = radii < .5
        margin = np.max(np.where(valid, mass-cost, -np.inf), axis=1)
        raw_margin = np.max(np.where(valid, mass-4*radii**2, -np.inf), axis=1)
        hit += int(np.count_nonzero(margin > 1e-12))
        internal += int(np.count_nonzero(raw_margin > 1e-12))
        near += int(np.count_nonzero(np.abs(margin) <= 1e-10))
        global_margin = max(global_margin, float(np.max(margin)))
    cell_area = ((w+1)/resolution)**2
    return {'n': n, 'resolution': resolution, 'sample_points': len(points),
            'halo_positive_samples': hit, 'internal_positive_samples': internal,
            'numerical_midpoint_halo_area': hit*cell_area,
            'numerical_midpoint_internal_area': internal*cell_area,
            'near_zero_best_margin_samples': near, 'max_sampled_mass_margin': global_margin,
            'pid': os.getpid(), 'elapsed_seconds': time.monotonic()-started}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--n', type=int, nargs='+', default=[5, 17, 65])
    parser.add_argument('--resolution', type=int, default=128)
    parser.add_argument('--workers', type=int, default=3)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    assert all(n >= 2 for n in args.n)
    started = time.monotonic()
    with concurrent.futures.ProcessPoolExecutor(max_workers=args.workers) as pool:
        results = list(pool.map(probe, [(n, args.resolution) for n in args.n]))
    out = {'schema': 'r072-diffuse-midpoint-probe-v1', 'pid': os.getpid(),
           'python': sys.version, 'numpy': np.__version__, 'platform': platform.platform(),
           'script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
           'seed': None, 'requested_workers': args.workers,
           'actual_returning_worker_pids': sorted({r['pid'] for r in results}),
           'pool_closed': True, 'results': results,
           'scope': 'Floating point samples only; midpoint sums have no certified error bound. '
                    'Cannot establish deletion area or finite-atomic residual transfer.',
           'count_rounding_tolerance': 1e-11, 'strict_margin_tolerance': 1e-12,
           'elapsed_seconds': time.monotonic()-started}
    args.output.write_text(json.dumps(out, indent=2)+'\n')
    print(json.dumps(out, indent=2))


if __name__ == '__main__':
    main()
