#!/usr/bin/env python3
"""Random-second-birth-rank importance design and guarded GPU dispatcher.

This file deliberately has no dependency on the live GPU source.  It derives
the *CLI domain* that makes the frozen r021 binary use a prescribed, unique
64-bit trajectory seed.  The layer-selection stream and trajectory-seed
stream have distinct domains.
"""
from fractions import Fraction
import argparse, hashlib, json, pathlib, subprocess, sys, tempfile

MASK = (1 << 64) - 1
GAMMA = 0x9E3779B97F4A7C15
GPU_BASE = 0xA022000000000000

def splitmix_word(x: int) -> int:
    """The bijective SplitMix64 output map at input x (mod 2^64)."""
    z = (x + GAMMA) & MASK
    z = ((z ^ (z >> 30)) * 0xBF58476D1CE4E5B9) & MASK
    z = ((z ^ (z >> 27)) * 0x94D049BB133111EB) & MASK
    return (z ^ (z >> 31)) & MASK

def q_table(N: int, K: int = 6):
    """Return (tau,qnum,qden), exactly q=3g/4+uniform/4."""
    if N < 2:
        raise ValueError("N must be at least 2")
    L = N - 1
    K = min(K, N - 2)
    tail = N - K - 1
    D = (1 << K) * tail
    den = 4 * D * L
    rows = []
    for tau in range(2, N + 1):
        j = tau - 1
        gnum = (1 << (K - j)) * tail if j <= K else 1
        qnum = 3 * L * gnum + D
        rows.append((tau, qnum, den))
    if sum(qnum for _, qnum, _ in rows) != den:
        raise AssertionError("q numerators do not sum to denominator")
    return rows

def below_rejection(state: int, den: int):
    """Uniform [0,den), consuming SplitMix words; no modulo bias."""
    if not 1 <= den <= 1 << 64:
        raise ValueError("den outside UInt64 range")
    lim = ((1 << 64) // den) * den
    rejects = 0
    while True:
        word = splitmix_word(state)
        state = (state + GAMMA) & MASK
        if word < lim:
            return state, word % den, rejects
        rejects += 1

def choose(rows, state: int):
    den = rows[0][2]
    state, u, rejects = below_rejection(state, den)
    c = 0
    for tau, num, _ in rows:
        c += num
        if u < c:
            return state, tau, num, den, rejects, u
    raise AssertionError("CDF did not cover accepted word")

def desired_seed(trajectory_domain: int, globalrep: int) -> int:
    """A bijective per-globalrep seed, provided globalrep values are unique."""
    if not 0 <= globalrep <= MASK:
        raise ValueError("globalrep outside UInt64")
    return splitmix_word((trajectory_domain + globalrep) & MASK)

def cli_domain_for_seed(seed: int, n: int, tau: int, localrep: int = 0) -> int:
    """Invert r021 CLI: seed=BASE+domain+n*1e6+tau*10000+localrep."""
    return (seed - GPU_BASE - n * 1_000_000 - tau * 10_000 - localrep) & MASK

def gpu_seed_from_domain(domain: int, n: int, tau: int, localrep: int = 0) -> int:
    return (GPU_BASE + domain + n * 1_000_000 + tau * 10_000 + localrep) & MASK

def rejection_certificate(den: int, bits: int = 64, enumerate_limit: int = 1 << 20):
    """Exact quotient/remainder certificate; enumerate only small intervals."""
    if den <= 0 or den > (1 << bits):
        raise ValueError("invalid denominator")
    lim = ((1 << bits) // den) * den
    ans = {"bits": bits, "den": den, "accept": lim,
           "reject": (1 << bits) - lim, "quotient": lim // den,
           "uniform_by_partition": True, "enumerated": lim <= enumerate_limit}
    if ans["enumerated"]:
        counts = [0] * den
        for x in range(lim): counts[x % den] += 1
        ans["counts"] = counts
        ans["uniform_by_enumeration"] = len(set(counts)) == 1
    return ans

def exact_n3(path: str):
    """Rational n=3 check: sum identity and random-tau variance identities."""
    x = json.load(open(path))
    counts = {int(k): int(v) for k, v in x['second_birth_rank_counts'].items()}
    sq = {int(k): int(v) for k, v in x['second_birth_rank_squared_totals'].items()}
    N = 9
    rows = q_table(N)
    A = sum(Fraction(counts.get(t, 0), 1) for t, _, _ in rows)
    DD = sum(Fraction(sq.get(t, 0), 1) for t, _, _ in rows)
    ez2 = sum(Fraction(counts.get(t, 0) ** 2 * den, num) for t, num, den in rows)
    ed2 = sum(Fraction(sq.get(t, 0) ** 2 * den, num) for t, num, den in rows)
    # Exact expectation under tau~q, now tests z_tau/q_tau rather than q=1.
    ez = sum(Fraction(aa, bb) * Fraction(counts.get(t, 0) * bb, aa)
             for t, aa, bb in rows)
    ed = sum(Fraction(aa, bb) * Fraction(sq.get(t, 0) * bb, aa)
             for t, aa, bb in rows)
    if ez != A or ed != DD:
        raise AssertionError("importance sum identity failed")
    return {
        "A2": str(A), "euclidean_squared_total": str(DD),
        "Q_euclidean": str(DD / (N * A)),
        "E_Z_over_q": str(ez), "E_D_over_q": str(ed),
        "Var_Z_over_q": str(ez2 - ez * ez),
        "Var_D_over_q": str(ed2 - ed * ed),
        "q": [{"tau": t, "qnum": aa, "qden": b} for t, aa, b in rows],
    }

def dispatch_records(n, reps, K, tau_domain, trajectory_domain):
    rows = q_table(n * n, K)
    state = tau_domain
    out, seen_seed = [], set()
    for globalrep in range(reps):
        state, tau, qnum, qden, rejects, u = choose(rows, state)
        seed = desired_seed(trajectory_domain, globalrep)
        cli_domain = cli_domain_for_seed(seed, n, tau)
        if gpu_seed_from_domain(cli_domain, n, tau) != seed or seed in seen_seed:
            raise RuntimeError("seed inversion/non-uniqueness failed")
        seen_seed.add(seed)
        out.append({"kind": "dispatch", "n": n, "globalrep": globalrep,
                    "tau": tau, "qnum": qnum, "qden": qden,
                    "tau_u": u, "tau_rejections": rejects,
                    "tau_domain_hex": f"{tau_domain:016x}",
                    "trajectory_domain_hex": f"{trajectory_domain:016x}",
                    "desired_seed_hex": f"{seed:016x}",
                    "domain_hex": f"{cli_domain:016x}", "localrep": 0})
    return rows, out

def run_gpu(binary, out_jsonl, n, P, moveL, method, rec):
    """One R=1 call, then insist that frozen CLI reports desired seed exactly."""
    with tempfile.TemporaryDirectory(dir=str(pathlib.Path(out_jsonl).parent)) as td:
        raw = pathlib.Path(td) / "one.jsonl"
        cmd = [binary, str(n), str(P), "1", str(rec['tau']), str(rec['tau']),
               str(raw), rec['domain_hex'], method, str(moveL)]
        cp = subprocess.run(cmd, check=False, text=True, capture_output=True)
        if cp.returncode:
            raise RuntimeError(f"GPU failure: {cp.returncode}: {cp.stderr[-1000:]}")
        got = [json.loads(s) for s in raw.read_text().splitlines() if s]
    reps = [r for r in got if r.get("kind") == "replicate"]
    if len(reps) != 1 or reps[0].get("seed_hex", "").lower().zfill(16) != rec['desired_seed_hex']:
        raise RuntimeError("frozen CLI seed did not match inverted desired seed")
    row = dict(rec); row.update(reps[0])
    row.update({"n": n, "P": P, "moveL": moveL, "method": method})
    row['gpu_localrep'] = row.pop('replicate')
    with open(out_jsonl, "a") as f:
        f.write(json.dumps(row, separators=(",", ":")) + "\n")

def main():
    p = argparse.ArgumentParser()
    p.add_argument('--N', type=int, help='number of graph vertices for design-only mode')
    p.add_argument('--n', type=int, help='side length for dispatch/GPU mode')
    p.add_argument('--K', type=int, default=6); p.add_argument('--reps', type=int, default=8)
    p.add_argument('--tau-domain', type=lambda s: int(s, 16), default=0x5210000000000000)
    p.add_argument('--trajectory-domain', type=lambda s: int(s, 16), default=0xA721000000000000)
    p.add_argument('--n3-exact'); p.add_argument('--gpu-binary'); p.add_argument('--out-jsonl')
    p.add_argument('--P', type=int, default=17); p.add_argument('--moveL', type=int, default=4)
    p.add_argument('--method', choices=('strat', 'multi'), default='strat')
    p.add_argument('--certificate-bits', type=int, default=16)
    a = p.parse_args()
    if (a.N is None) == (a.n is None): p.error('give exactly one of --N or --n')
    N = a.N if a.N is not None else a.n * a.n
    rows = q_table(N, a.K)
    # This certificate is O(1), so n=24/32 never tries an immense counts
    # array. The separate small test below exercises actual rejection.
    cert = rejection_certificate(rows[0][2], 64)
    # A non-divisor test exercises the rejection branch even when the chosen
    # small-n q denominator happens to divide the simulated word range.
    adversarial_cert = rejection_certificate(1001, a.certificate_bits)
    if not cert['uniform_by_partition'] or not adversarial_cert.get('uniform_by_enumeration'):
        raise RuntimeError('rejection certificate failed')
    ans = {"kind": "random_tau_design", "N": N, "K": min(a.K, N - 2),
           "qden": rows[0][2], "rejection_certificate": cert,
           "adversarial_rejection_certificate": adversarial_cert}
    if a.n3_exact: ans['n3_exact'] = exact_n3(a.n3_exact)
    if a.n is not None:
        _, recs = dispatch_records(a.n, a.reps, a.K, a.tau_domain, a.trajectory_domain)
        ans['dispatch'] = recs
        if a.gpu_binary or a.out_jsonl:
            if not (a.gpu_binary and a.out_jsonl): p.error('--gpu-binary and --out-jsonl together')
            out = pathlib.Path(a.out_jsonl)
            if out.exists(): raise RuntimeError(f'refusing to append to existing {out}')
            binary_path = pathlib.Path(a.gpu_binary).resolve()
            manifest = dict(ans); manifest.update({"n": a.n, "reps": a.reps, "P": a.P, "moveL": a.moveL,
                "method": a.method, "gpu_binary": str(binary_path),
                "gpu_binary_sha256": hashlib.sha256(binary_path.read_bytes()).hexdigest(),
                "dispatcher_sha256": hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest(),
                "python": sys.version, "invocation": sys.argv})
            out.with_suffix(out.suffix + '.manifest.json').write_text(json.dumps(manifest, indent=2) + "\n")
            for rec in recs: run_gpu(a.gpu_binary, str(out), a.n, a.P, a.moveL, a.method, rec)
    print(json.dumps(ans, indent=2))

if __name__ == '__main__': main()
