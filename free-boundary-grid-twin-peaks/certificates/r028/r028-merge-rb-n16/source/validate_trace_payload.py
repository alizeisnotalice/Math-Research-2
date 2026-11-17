#!/usr/bin/env python3
"""Structural P17 trace-payload audit for the frozen r026 derivative.

This script independently recomputes the RB payload aggregates from the
trace's terminal payload.  It deliberately validates finite trace fields only,
not mixing, SMC accuracy, or an asymptotic cone statement.
"""
import hashlib
import json
import math
import sys
from pathlib import Path

SOURCE = '3a9f3c0a8c10409346c25f16f2d82463f928bfaacf7307b4a88a770ec8f69ec6'
BASE = 0xA022000000000000

def main():
    raw, trace, source = map(Path, sys.argv[1:4])
    if hashlib.sha256(source.read_bytes()).hexdigest() != SOURCE:
        raise SystemExit('calibration derivative source hash mismatch')
    lines = [json.loads(line) for line in raw.read_text().splitlines()]
    if len(lines) != 4 or lines[0].get('kind') != 'metadata' or lines[-1].get('kind') != 'run_end':
        raise SystemExit('expected metadata, replicate, RB payload, run_end')
    meta, rep, rb, end = lines
    if meta.get('implementation') != 'r026-tau2-firstmerge-rb-terminal':
        raise SystemExit('unexpected raw metadata')
    n, P, R, lo, hi, L, method = (meta.get(k) for k in ('n', 'particles', 'replicates', 'tau_lo', 'tau_hi', 'move_steps', 'method'))
    if not (isinstance(n, int) and 3 <= n <= 32 and (P, R, lo, hi, L, method) == (17, 1, 2, 2, 16, 'strat')):
        raise SystemExit('unexpected trace metadata')
    if rep.get('kind') != 'replicate' or rb.get('kind') != 'firstmerge_rb_calibration' or rep.get('replicate') != rb.get('replicate') != 0 or rep.get('dead'):
        raise SystemExit('replicate/RB pairing failure')
    domain = int(meta['domain'], 16)
    seed = (BASE + domain + n * 1_000_000 + 2 * 10_000) & ((1 << 64) - 1)
    if int(rep['seed_hex'], 16) != seed:
        raise SystemExit('seed mismatch')
    data = json.loads(trace.read_text())
    if (data.get('implementation'), data.get('n'), data.get('P'), data.get('tau'), data.get('L'), data.get('method'), int(data.get('seed_hex'), 16)) != ('r026-tau2-firstmerge-rb-terminal', n, P, 2, L, method, seed):
        raise SystemExit('trace metadata mismatch')
    terminal = [row for row in data.get('rows', []) if row.get('terminal')]
    if len(terminal) != 1 or terminal[0].get('first_merge_rb_bad') != 0:
        raise SystemExit('terminal payload guard failure')
    payload = terminal[0].get('first_merge_rb_payload_ab_r1_r2')
    if not isinstance(payload, list) or len(payload) != 4 * P:
        raise SystemExit('payload length mismatch')
    out = [0.0] * 4
    for i in range(P):
        a, b, x, y = payload[4*i:4*i+4]
        if not (isinstance(a, int) and isinstance(b, int) and a > 0 and b > 0 and 0 <= x < n*n and 0 <= y < n*n):
            raise SystemExit('invalid payload item')
        q = a + b
        m = q * (q - 1) / (2 * a * b)
        dx, dy = x % n - y % n, x // n - y // n
        d2 = dx * dx + dy * dy
        ax, ay, bx, by = 2*(x % n)-(n-1), 2*(x//n)-(n-1), 2*(y % n)-(n-1), 2*(y//n)-(n-1)
        f, g = ax*ax + ay*ay + bx*bx + by*by, ax*bx + ay*by
        if f - 2*g != 4*d2:
            raise SystemExit('pointwise root identity failure')
        for j, value in enumerate((m, m*d2, m*f, m*g)):
            out[j] += value / P
    for key, value in zip(('mean_m', 'mean_m_d2', 'mean_m_F0', 'mean_m_G0'), out):
        if not math.isfinite(rb.get(key, float('nan'))) or abs(rb[key] - value) > 1e-12:
            raise SystemExit(f'payload mean mismatch: {key}')
    if abs(rb['mean_m_F0'] - 2*rb['mean_m_G0'] - 4*rb['mean_m_d2']) > 1e-10:
        raise SystemExit('RB aggregate identity failure')
    print(json.dumps({'audit_passed': True, 'n': n, 'P': P, 'scope': 'P17 trace and RB payload structural audit only; not an accuracy, mixing, or cone claim'}))

if __name__ == '__main__':
    main()
