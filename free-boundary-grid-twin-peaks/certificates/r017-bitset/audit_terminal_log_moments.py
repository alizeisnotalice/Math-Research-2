#!/usr/bin/env python3
"""Audit terminal serialization and UInt64 bounds independently of Swift."""
import json
import math
from pathlib import Path

ROOT = Path(__file__).parent


def replicate(name):
    rows = [json.loads(x) for x in (ROOT / name).read_text().splitlines()]
    return next(row for row in rows if row.get("kind") == "replicate")


def logsumexp(log_weights):
    """Required aggregation for non-dead tau terms; never sum underflowed Z."""
    finite = [x for x in log_weights if math.isfinite(x)]
    pivot = max(finite)
    return pivot + math.log(sum(math.exp(x - pivot) for x in finite))


n12 = replicate("frozen-final-n12-tau2-p17-domain1.jsonl")
n32_mid = replicate("frozen-final-n32-tau512-p17-domain3.jsonl")
n32_nm1 = replicate("frozen-final-n32-tau1023-p17-domain3.jsonl")
n32_n = replicate("frozen-final-n32-tau1024-p17-domain3.jsonl")
n32_one = replicate("frozen-final-n32-onepeak-p17-domain20.jsonl")

# Positive-range identity: direct exp(logZ) and emitted Z agree when representable.
assert not n12["dead"] and math.isfinite(n12["logZ"])
assert n12["Z"] == math.exp(n12["logZ"])
assert n12["U"] == n12["Z"] * n12["U_over_Z"]

# Death has null logZ, whereas an underflowed but live run retains finite logZ.
for row in (n32_mid, n32_one):
    assert not row["dead"] and math.isfinite(row["logZ"]) and row["Z"] == 0.0
for row in (n32_nm1, n32_n):
    assert row["dead"] and row["logZ"] is None and row["Z"] == 0.0
assert n32_one["U"] is None and n32_one["U_over_Z"] is None and n32_one["scaled_U"] is None

# Integer terminal accumulator: each particle has at most N terms and
# squared Manhattan-coordinate displacement at most 2(n-1)^2.
# The group kernel sums at most 256 particles.  The P^2*N guard gives
# P <= floor(sqrt((2^64-1)/N)); hence the total raw sum is safely UInt64.
N = 32 * 32
P_max = math.isqrt(((1 << 64) - 1) // N)
per_particle = N * 2 * (32 - 1) ** 2
per_group = 256 * per_particle
raw_total_bound = P_max * per_particle
assert per_group < 1 << 64 and raw_total_bound < 1 << 64

# Example of the only valid cross-tau normalizer operation.  These values are
# intentionally not a probability estimate; they demonstrate log-domain form.
demo_log_normalizer = logsumexp([n32_mid["logZ"], n32_one["logZ"]])

result = {
    "status": "passed",
    "n12_positive_range_identity": {k: n12[k] for k in ("Z", "U", "logZ", "U_over_Z", "dead")},
    "n32_live_underflow": {
        "tau512": {k: n32_mid[k] for k in ("Z", "logZ", "dead")},
        "tau1025_onepeak": {k: n32_one[k] for k in ("Z", "U", "logZ", "U_over_Z", "scaled_U", "dead")},
    },
    "n32_dead": {
        "tau1023": {k: n32_nm1[k] for k in ("Z", "logZ", "dead")},
        "tau1024": {k: n32_n[k] for k in ("Z", "logZ", "dead")},
    },
    "integer_bounds": {"N": N, "guard_max_P": P_max, "per_particle": per_particle,
                       "per_group": per_group, "raw_total_bound": raw_total_bound},
    "logsumexp_demo": demo_log_normalizer,
    "aggregation_rule": "For aligned non-dead tau terms, use ell=logsumexp(logZ_tau) and Q=sum(exp(logZ_tau-ell)*q_tau); never sum emitted Z values.",
}
(ROOT / "terminal-log-moment-audit.json").write_text(json.dumps(result, indent=2) + "\n")
print(json.dumps(result, indent=2))
