"""Run exact CRT on selected compact regions with a bounded job queue."""
import argparse
from datetime import datetime, timezone
import hashlib
import json
from math import comb, factorial, gcd, prod
from pathlib import Path
import re
import signal
import subprocess
import os
import sys
import time


HERE = Path(__file__).resolve().parent
MIB = 1 << 20


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verified_primes(path):
    certificate = json.loads(path.read_text())
    known = set()

    def visit(p):
        if p in known:
            return
        node = certificate['nodes'][str(p)]
        assert int(node['p']) == p
        if p == 2:
            assert node == {'p': '2', 'base': 'two'}
        else:
            factors = {int(q): k for q, k in node['p_minus_1_factorization'].items()}
            assert all(2 <= q < p and k >= 1 for q, k in factors.items())
            assert prod(q**k for q, k in factors.items()) == p - 1
            for q in factors:
                visit(q)
            a = int(node['witness'])
            assert 1 < a < p and pow(a, p - 1, p) == 1
            assert all(gcd(pow(a, (p - 1) // q, p) - 1, p) == 1 for q in factors)
        known.add(p)

    primes = tuple(map(int, certificate['roots']))
    assert len(primes) == len(set(primes)) == 3
    for p in primes:
        assert p.bit_length() == 61
        visit(p)
    assert prod(primes) > factorial(41)
    return primes, len(known)


def independent_sizes(profile):
    """Forward cardinality DP, independent of the C++ suffix DP."""
    last = {0: [1]}
    for width in map(int, profile):
        after = {}
        for s in range(1 << width):
            if s & (s << 1):
                continue
            shift = bin(s).count('1')
            current = [0] * (max(map(len, last.values())) + shift)
            for t, counts in last.items():
                if s & t:
                    continue
                for k, count in enumerate(counts):
                    current[k + shift] += count
            after[s] = current
        last = after
    histogram = [0] * max(map(len, last.values()))
    for counts in last.values():
        for k, count in enumerate(counts):
            histogram[k] += count
    while not histogram[-1]:
        histogram.pop()
    return histogram


def safe_roots(profile):
    n = len(profile)
    region = {y * n + x for y, w in enumerate(map(int, profile)) for x in range(w)}
    out = set()
    for v in region:
        x, y = v % n, v // n
        if x >= (n + 3) // 4:
            continue
        closed = {Y * n + X for X, Y in ((x - 1, y), (x + 1, y), (x, y - 1), (x, y + 1))
                  if 0 <= X < n and 0 <= Y < n}
        if closed <= region:
            out.add(v)
    return out


def crt(residues, primes):
    result, modulus = 0, 1
    for r, p in zip(residues, primes):
        result += ((r - result) * pow(modulus, -1, p) % p) * modulus
        modulus *= p
    assert all(result % p == r for r, p in zip(residues, primes))
    return result


def command(h, ix, budget):
    return [str(HERE / 'root_mod'), h, str(ix), '--max-table-mib', str(budget)]


def plan(h, budget):
    p = subprocess.run(command(h, 0, budget) + ['--plan'], capture_output=True, text=True, timeout=5, check=True)
    record = json.loads(p.stdout)
    assert p.stderr == '' and record['mode'] == 'plan' and record['profile'] == h
    assert record['states_by_layer'] == independent_sizes(h)
    assert record['independent_sets'] == sum(record['states_by_layer'])
    assert record['allowed_roots'] == len(safe_roots(h))
    assert record['area'] == sum(map(int, h))
    return record


def stop_group(proc):
    if proc.poll() is not None:
        return
    try:
        os.killpg(proc.pid, signal.SIGTERM)
    except ProcessLookupError:
        pass
    try:
        proc.wait(timeout=3)
    except subprocess.TimeoutExpired:
        try:
            os.killpg(proc.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
        proc.wait()


def run_queue(output, planned, primes, workers, table_budget, aggregate_budget, timeout):
    pending = [(h, ix) for h in planned for ix in range(len(primes))]
    active = {}
    done = {}
    max_active_payload = 0
    completed = output / 'receipts.jsonl'
    with completed.open('w') as receipt:
        try:
            while pending or active:
                used = sum(planned[h]['max_table_bytes'] for h, _ in active)
                while pending and len(active) < workers:
                    h, ix = pending[0]
                    payload = planned[h]['max_table_bytes']
                    if payload + used > aggregate_budget:
                        break
                    pending.pop(0)
                    stem = output / f'{h}-p{ix}'
                    stdout, resource = stem.with_suffix('.json'), stem.with_suffix('.resource')
                    with stdout.open('w') as out, resource.open('w') as err:
                        argv = ['/usr/bin/time', '-l', *command(h, ix, table_budget), '--progress']
                        proc = subprocess.Popen(argv, stdout=out, stderr=err, start_new_session=True)
                    active[h, ix] = (proc, stdout, resource, time.monotonic(), argv)
                    used += payload
                    max_active_payload = max(max_active_payload, used)
                if pending and not active:
                    raise ValueError('pending job cannot fit aggregate table budget')
                for key, active_job in list(active.items()):
                    proc, stdout, resource, started, argv = active_job
                    status = proc.poll()
                    if status is None:
                        if time.monotonic() - started > timeout:
                            raise TimeoutError(f'timing out job {key}')
                        continue
                    del active[key]
                    if status != 0:
                        raise RuntimeError(f'job {key} exited with {status}; see {resource}')
                    h, ix = key
                    data = json.loads(stdout.read_text())
                    assert data['mode'] == 'count' and data['profile'] == h
                    assert data['prime_index'] == ix and int(data['prime']) == primes[ix]
                    assert data['processed_by_layer'] == data['states_by_layer'] == planned[h]['states_by_layer']
                    assert data['max_table_bytes'] == planned[h]['max_table_bytes']
                    assert len(data['root_residues']) == len(h)**2
                    assert all(0 <= int(x) < primes[ix] for x in data['root_residues'])
                    timing = resource.read_text()
                    rss = re.search(r'^\s*(\d+)\s+maximum resident set size$', timing, re.M)
                    assert rss is not None, timing
                    record = {'profile': h, 'prime_index': ix, 'pid': proc.pid, 'exit_code': status,
                              'max_rss_bytes': int(rss[1]), 'wall_seconds': time.monotonic() - started,
                              'stdout_sha256': digest(stdout), 'resource_sha256': digest(resource), 'argv': argv}
                    print(json.dumps(record), file=receipt, flush=True)
                    done[key] = data
                    print('completed', h, ix, 'RSS bytes', rss[1], flush=True)
                if active:
                    time.sleep(0.15)
        finally:
            for proc, _, _, _, _ in active.values():
                stop_group(proc)
    return done, max_active_payload


def build_outputs(profiles, raw, primes):
    exact = {}
    totals = {}
    for h in profiles:
        vectors = [list(map(int, raw[h, ix]['root_residues'])) for ix in range(len(primes))]
        counts = [crt(rs, primes) for rs in zip(*vectors)]
        assert {v for v, count in enumerate(counts) if count > 0} == safe_roots(h)
        assert sum(counts) <= factorial(sum(map(int, h)))
        exact[h] = counts
        totals[h] = sum(counts)
    words = {}
    for h in profiles:
        mirror = ''.join(str(len(h)-int(c)) for c in h)
        if mirror in totals:
            words[h] = str(comb(len(h)**2, sum(map(int, h))) * totals[h] * totals[mirror])
    return {h: list(map(str, counts)) for h, counts in exact.items()}, {h: str(v) for h, v in totals.items()}, words


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('output', type=Path)
    parser.add_argument('profiles', nargs='+')
    parser.add_argument('--workers', type=int, default=3, choices=range(1, 4))
    parser.add_argument('--max-table-mib', type=int, default=1024)
    parser.add_argument('--aggregate-table-mib', type=int, default=3072)
    parser.add_argument('--timeout', type=int, default=1800)
    parser.add_argument('--plan-only', action='store_true')
    options = parser.parse_args()
    assert len(options.profiles) == len(set(options.profiles))
    assert min(options.max_table_mib, options.aggregate_table_mib, options.timeout) > 0
    root_cert = HERE / 'pratt-61bit.json'
    primes, nodes = verified_primes(root_cert)
    planned = {h: plan(h, options.max_table_mib) for h in options.profiles}
    options.output.mkdir(exist_ok=False)
    metadata = {'schema': 'a022-root-crt-batch-v1', 'python': sys.version,
                'started_utc': datetime.now(timezone.utc).isoformat(), 'primes': list(map(str, primes)),
                'pratt_verified_nodes': nodes, 'prime_cert_sha256': digest(root_cert),
                'counter_source_sha256': digest(HERE / 'root_mod.cpp'), 'counter_binary_sha256': digest(HERE / 'root_mod'),
                'driver_sha256': digest(Path(__file__)), 'workers': options.workers,
                'aggregate_table_budget_mib': options.aggregate_table_mib, 'per_job_table_budget_mib': options.max_table_mib,
                'planned': planned}
    (options.output / 'plan.json').write_text(json.dumps(metadata, indent=2) + '\n')
    if options.plan_only:
        print('planned', len(planned), 'profiles')
        return
    assert all(p['fits_budget'] for p in planned.values())
    raw, simultaneous = run_queue(options.output, planned, primes, options.workers, options.max_table_mib,
                                  options.aggregate_table_mib * MIB, options.timeout)
    roots, totals, events = build_outputs(options.profiles, raw, primes)
    summary = dict(metadata, scope='Rooted one-peak regional counts and individual cut counts. No A2 normalization.',
                   completed_jobs=len(raw), max_active_table_payload_bytes=simultaneous,
                   root_counts=roots, safe_root_totals=totals, individual_cut_counts=events,
                   finished_utc=datetime.now(timezone.utc).isoformat())
    (options.output / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')


if __name__ == '__main__':
    def on_term(_signum, _frame):
        raise KeyboardInterrupt('terminated')
    signal.signal(signal.SIGTERM, on_term)
    main()
