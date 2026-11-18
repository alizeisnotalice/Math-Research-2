"""Positive local-density cuts for the symmetric quartic contact."""
import argparse
import hashlib
import json
import math
import os
import platform
import resource
import socket
import time
from pathlib import Path

import numpy as np
import scipy
from scipy.special import logsumexp, roots_legendre
from scipy.stats import qmc

PI2 = 2 * np.pi
SCALES = np.array([0.03125, 0.125, 0.5, 1., 2., 8., 32.])
SHOULDERS = np.array([-12., -5., -2., -1., 0., 1., 2., 5., 12.])


def log_integrals(c, x, m, eps, delta, order, radius):
    """Compute normalized masses for cores cut at radius and at v=M."""
    nodes, weights = roots_legendre(order)
    scale_h = np.sqrt(eps)
    log_target = -0.5*(np.log(PI2*(m+eps)) + (x+delta)**2/(m+eps))
    upto_radius = []
    upto_cap = []
    for sign, k, power in [(-1., c, 4), (1., c, 4)]:
        cap = np.nextafter((m/k)**(1./power), 0.)
        cap = np.minimum(cap, (m/k)**(1./power)*(1.-8*np.finfo(float).eps))
        assert np.all(k*cap**power <= m)
        small = np.minimum(cap, radius)
        shoulder = np.maximum(0., sign*delta)
        width = np.sqrt(eps+k*shoulder**power)
        exponent_one = 1./np.sqrt(k)
        groups = [np.zeros_like(cap[:, None]), small[:, None], cap[:, None],
                  scale_h[:, None]*SCALES,
                  (eps/k)[:, None]**(1./power)*SCALES,
                  exponent_one[:, None]*SCALES,
                  ((delta*delta/k)**.25)[:, None] * SCALES,
                  shoulder[:, None]+width[:, None]*SHOULDERS,
                  (-sign*x)[:, None]+SHOULDERS]
        edges = np.sort(np.clip(np.concatenate(groups, axis=1), 0., cap[:, None]), axis=1)
        midpoint = (edges[:, 1:] + edges[:, :-1])*0.5
        half_width = (edges[:, 1:] - edges[:, :-1])*0.5
        u = midpoint[:, :, None] + half_width[:, :, None]*nodes
        var = eps[:, None, None] + k[:, None, None]*u**power
        log_kernel = (-np.log(PI2) - 0.5*(np.log(var)+(x[:, None, None]+sign*u)**2
                                           +(delta[:, None, None]-sign*u)**2/var)
                      - log_target[:, None, None])
        with np.errstate(divide="ignore"):
            log_term = log_kernel+np.log(half_width[:, :, None])+np.log(weights)
        panel = logsumexp(log_term, axis=2)
        upto_cap.append(logsumexp(panel, axis=1))
        upto_radius.append(logsumexp(np.where(edges[:, 1:] <= small[:, None], panel, -np.inf), axis=1))
    return (np.logaddexp(*upto_radius), np.logaddexp(*upto_cap))


def samples(shard, shards, exponent):
    grid = qmc.Sobol(5, scramble=True, rng=np.random.default_rng(821602+shard)).random_base2(exponent)
    a, b, e, y, g = grid.T
    p = (shard + a)/shards
    c = np.exp2(-10.+22.*p)
    log_m = b*np.log(12.)
    m = np.exp(log_m)
    threshold = np.sqrt((m*log_m)/np.expm1(log_m))
    x = threshold + np.exp2(-17.+19.*g)
    assert np.all(x > threshold)
    eps = np.exp2(-84.+88.*e)
    offset_mode = np.arange(len(a)) % 5
    delta = np.select([offset_mode == 0, offset_mode == 1, offset_mode == 2, offset_mode == 3],
                      [0., (y-.5)*8.*np.sqrt(eps), (y-.5)*20.*np.sqrt(eps), (y-.5)*.5],
                      default=(y-.5)*8.)
    return c, x, m, eps, delta


def take_top(ratios, m, eps, delta, c, x, n=12):
    indices = np.argpartition(ratios, -n)[-n:]
    indices = sorted(indices, key=lambda i: ratios[i], reverse=True)
    return [dict(index=int(i), c=float(c[i]), x=float(x[i]),
                 m=float(m[i]), eps=float(eps[i]), delta=float(delta[i]), log_ratio=float(ratios[i]))
            for i in indices]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("shard", type=int)
    parser.add_argument("--shards", type=int, default=15)
    parser.add_argument("--exponent", type=int, default=16)
    parser.add_argument("--order", type=int, default=24)
    parser.add_argument("--radius", type=float, default=.25)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    start = time.monotonic()
    c, x, m, eps, delta = samples(args.shard, args.shards, args.exponent)
    local = np.empty(len(c))
    cap = np.empty(len(c))
    for lo in range(0, len(c), 192):
        i = slice(lo, lo+192)
        local[i], cap[i] = log_integrals(c[i], x[i], m[i], eps[i], delta[i], args.order, args.radius)
    assert np.all(np.isfinite(local)) and np.all(np.isfinite(cap))
    assert np.all(local <= cap + 1e-12)
    summary = {}
    for name, values in [("radius_capped", local), ("variance_capped", cap)]:
        top = take_top(values, m, eps, delta, c, x)
        i = np.array([r["index"] for r in top])
        for resolution in (48, 96):
            local_ref, cap_ref = log_integrals(c[i], x[i], m[i], eps[i], delta[i],
                                               resolution, args.radius)
            ref = local_ref if name == "radius_capped" else cap_ref
            for result, value in zip(top, ref, strict=True):
                result[f"log_ratio_order_{resolution}"] = float(value)
        summary[name] = dict(positive_cuts=int(np.count_nonzero(values > 0)), top=top)
    report = dict(shard=args.shard, shards=args.shards, exponent=args.exponent, points=len(c),
                  range_c=[min(c), max(c)], range_x=[min(x), max(x)], order=args.order, fixed_radius=args.radius,
                  host=socket.gethostname(), machine=platform.machine(),
                  python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__,
                  process_id=os.getpid(), elapsed_seconds=time.monotonic()-start,
                  rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                  cpu_seconds=time.process_time(), script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  numerical_summary=summary)
    args.out.write_text(json.dumps(report, indent=2)+"\n")
    print(json.dumps(dict(shard=args.shard, points=len(c), secs=report["elapsed_seconds"], pid=os.getpid())))


if __name__ == "__main__":
    main()
