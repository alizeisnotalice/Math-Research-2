"""Plan or count top-two protected peak bounds with four certified moduli."""
import argparse
from concurrent.futures import ThreadPoolExecutor
from decimal import Decimal, localcontext
from fractions import Fraction
import hashlib
import json
from math import comb, factorial, gcd, prod
import os
from pathlib import Path
import re
import signal
import subprocess
import sys

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
    assert len(primes) == len(set(primes)) == 4
    for p in primes:
        assert p.bit_length() == 61
        visit(p)
    assert prod(primes) > factorial(50)
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


def crt(residues, primes):
    result, modulus = 0, 1
    for r, p in zip(residues, primes):
        result += ((r - result) * pow(modulus, -1, p) % p) * modulus
        modulus *= p
    assert all(result % p == r for r, p in zip(residues, primes))
    return result


def complement(h):
    return ''.join(str(len(h)-int(c)) for c in h)


def command(profile, ix, budget):
    return [str(HERE / 'protected_mod'), profile, str(ix), '--max-table-mib', str(budget)]


def plan(profile, budget, primes):
    done = subprocess.run(command(profile, 0, budget) + ['--plan'], stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE, text=True, check=True, timeout=5)
    out = json.loads(done.stdout)
    assert done.stderr == '' and out['mode'] == 'plan'
    assert out['profile'] == profile and out['area'] == sum(map(int, profile))
    assert int(out['prime']) == primes[0]
    assert out['protected_vertices'] == out['area'] - len(profile)
    sizes = independent_sizes(''.join(str(int(c)-1) for c in profile))
    assert sizes == out['states_by_layer'] and sum(sizes) == out['independent_sets']
    slots = [0] + [(0 if m == 0 else 1 << (((4*m+2)//3)-1).bit_length()) for m in sizes[1:]]
    assert out['slots_by_layer'] == slots
    peak = max(16 * (a+b) for a,b in zip(slots, slots[1:])) if len(slots) > 1 else 0
    assert out['max_table_bytes'] == peak
    assert out['fits_budget'] == (peak <= budget*MIB)
    return out


def stop_group(proc):
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


def count(profile, idx, prime, output, planned, budget, timeout):
    raw = output / f'{profile}-p{idx}.json'
    usage = output / f'{profile}-p{idx}.resource'
    cmd = ['/usr/bin/time', '-l', *command(profile, idx, budget), '--progress']
    with raw.open('x') as stdout, usage.open('x') as stderr:
        proc = subprocess.Popen(cmd, stdout=stdout, stderr=stderr, start_new_session=True)
        try:
            status = proc.wait(timeout=timeout)
        finally:
            if proc.poll() is None:
                stop_group(proc)
    if status:
        raise RuntimeError(f'{profile} prime {idx} exited {status}, see {usage}')
    data = json.loads(raw.read_text())
    assert data['mode'] == 'count' and data['profile'] == profile
    assert int(data['prime']) == prime and data['prime_index'] == idx
    assert data['processed_by_layer'] == data['states_by_layer'] == planned['states_by_layer']
    assert data['max_table_bytes'] == planned['max_table_bytes']
    fields = ('peak_frequency_residues', 'unprotected_max_frequency_residues', 'bad_second_frequency_residues')
    for field in fields:
        assert len(data[field]) == len(planned['states_by_layer'])
        assert all(0 <= x < prime for x in data[field])
    rss = re.search(r'^\s*(\d+)\s+maximum resident set size$', usage.read_text(), re.M)
    assert rss is not None
    receipt = {'profile': profile, 'prime_index': idx, 'pid': proc.pid, 'exit_code': status,
               'max_rss_bytes': int(rss[1]), 'raw_sha256': digest(raw), 'resource_sha256': digest(usage), 'argv': cmd}
    return (profile, idx), data, receipt


def frequencies(profile, parts, primes):
    fields = ('peak_frequency_residues', 'unprotected_max_frequency_residues', 'bad_second_frequency_residues')
    f, t, v = [[crt(values, primes) for values in zip(*(parts[profile, i][s] for i in range(len(primes))))]
               for s in fields]
    a = sum(map(int, profile))
    assert sum(f) == factorial(a) and sum(t) == len(profile)*factorial(a-1)
    assert all(0 <= y <= x and 0 <= z <= x-y for x,y,z in zip(f,t,v))
    return f, t, v


def merge(profile, regional):
    f, t, v = [a + [0]*max(0, 3-len(a)) for a in regional[profile]]
    g, u, w = [a + [0]*max(0, 3-len(a)) for a in regional[complement(profile)]]
    a = sum(map(int, profile))
    total = len(profile)**2
    all_z_le2 = comb(total, a)*sum(f[i]*g[j] for i in range(3) for j in range(3-i))
    exclude_first = sum(comb(total-1, a-1)*t[i]*g[2-i] + comb(total-1, a)*f[i]*u[2-i] for i in range(3))
    exclude_second = sum(comb(total-2, a-2)*v[i]*g[2-i] + comb(total-2, a)*f[i]*w[2-i] +
                         comb(total-2, a-1)*((f[i]-t[i])*u[2-i] + t[i]*(g[2-i]-u[2-i])) for i in range(3))
    result = all_z_le2-exclude_first-exclude_second
    assert 0 <= result <= all_z_le2 <= factorial(total)
    ratio = Fraction(result, factorial(total))
    with localcontext() as ctx:
        ctx.prec = 55
        decimal = str(Decimal(ratio.numerator)/Decimal(ratio.denominator))
    return {'cut': profile, 'reflected_other_side': complement(profile), 'safe_peak_at_most_two_count': str(all_z_le2),
            'global_max_excluded_count': str(exclude_first), 'global_second_excluded_count': str(exclude_second),
            'A2_upper_count': str(result), 'A2_upper_over_factorial_exact': str(ratio), 'A2_upper_over_factorial_decimal': decimal}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('output', type=Path)
    ap.add_argument('cuts', nargs='*', default=['5555555555'])
    ap.add_argument('--max-table-mib', type=int, default=512)
    ap.add_argument('--aggregate-table-mib', type=int, default=2048)
    ap.add_argument('--workers', type=int, choices=range(1, 5), default=2)
    ap.add_argument('--timeout', type=int, default=900)
    ap.add_argument('--execute', action='store_true', help='run enumeration; default produces only the resource plan')
    cli = ap.parse_args()
    assert len(set(cli.cuts)) == len(cli.cuts) and cli.cuts
    assert min(cli.timeout, cli.max_table_mib, cli.aggregate_table_mib) > 0
    primes, certificate_nodes = verified_primes(HERE / 'pratt-61bit.json')
    profiles = list(dict.fromkeys(p for h in cli.cuts for p in (h, complement(h))))
    planned = {p: plan(p, cli.max_table_mib, primes) for p in profiles}
    largest = max(q['max_table_bytes'] for q in planned.values())
    cli.output.mkdir(exist_ok=False)
    meta = {'schema': 'a022-four-mod-protected-peak-bound-v1', 'scope': 'Exact finite A2 upper bounds',
            'prime_certificate_nodes': certificate_nodes, 'primes': list(map(str, primes)), 'python': sys.version,
            'prime_cert_sha256': digest(HERE/'pratt-61bit.json'), 'source_sha256': digest(HERE/'protected_mod.cpp'),
            'driver_sha256': digest(Path(__file__)), 'planned': planned, 'cuts': cli.cuts, 'workers': cli.workers,
            'max_concurrent_table_payload_bytes': largest*min(cli.workers, len(planned)*len(primes)),
            'regional_jobs': len(planned)*len(primes)}
    (cli.output/'plan.json').write_text(json.dumps(meta,indent=2)+'\n')
    if not cli.execute:
        print(f'plan-only, {meta["regional_jobs"]} jobs, largest worker table {largest/MIB:g} MiB')
        return
    assert all(q['fits_budget'] for q in planned.values())
    assert meta['max_concurrent_table_payload_bytes'] <= cli.aggregate_table_mib*MIB
    jobs = [(h,i,p) for h in profiles for i,p in enumerate(primes)]
    with ThreadPoolExecutor(max_workers=cli.workers) as pool:
        observed = list(pool.map(lambda job: count(*job, cli.output, planned[job[0]], cli.max_table_mib, cli.timeout), jobs))
    residues = {key: value for key,value,_ in observed}
    histograms = {h: frequencies(h, residues, primes) for h in profiles}
    exact = [merge(h, histograms) for h in cli.cuts]
    ledger = dict(meta, completed_jobs=len(observed), max_worker_rss_bytes=max(r['max_rss_bytes'] for _,_,r in observed),
                  histogram_field_order=['all', 'unsafe_max', 'unsafe_nonadjacent_second_after_safe_max'],
                  region_histograms={h: [list(map(str, row)) for row in stat] for h,stat in histograms.items()},
                  upper_counts=exact, receipts=[r for _,_,r in observed])
    (cli.output/'summary.json').write_text(json.dumps(ledger, indent=2)+'\n')
    for row in exact:
        print(row['cut'], row['A2_upper_count'], row['A2_upper_over_factorial_decimal'])


if __name__ == '__main__':
    main()
