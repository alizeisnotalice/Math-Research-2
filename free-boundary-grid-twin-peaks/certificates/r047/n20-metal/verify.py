"""Check a completed tau=2 Metal batch and aggregate by unnormalized mass."""
import hashlib
import json
import math
from pathlib import Path
import sys


def sha(file):
    return hashlib.sha256(file.read_bytes()).hexdigest()


def check(raw_path, receipt_path, resource_path, source_path):
    receipt = json.loads(receipt_path.read_text())
    assert receipt["raw_sha256"] == sha(raw_path)
    assert receipt["resource_sha256"] == sha(resource_path)
    assert receipt["source_sha256"] == sha(source_path)
    assert receipt["returncode"] == receipt["caffeinate_returncode"] == 0
    assert receipt["waited_reaped"] and not receipt["remote_process_presence_after_wait"]
    raw = [json.loads(line) for line in raw_path.read_text().splitlines()]
    m = raw[0]
    expected = receipt["parameters"]
    n, P, R = (m[k] for k in ("n", "particles", "replicates"))
    assert (n, P, R) == tuple(expected[k] for k in ("n", "P", "R"))
    assert m["kind"] == "metadata" and m["implementation"] == "r026-tau2-firstmerge-rb-terminal"
    assert m["tau_lo"] == m["tau_hi"] == expected["tau"] == 2
    assert m["move_steps"] == expected["L"] and m["domain"] == expected["domain"]
    assert m["method"] == "strat"
    assert raw[-1]["kind"] == "run_end" and len(raw) == 2 * R + 2 == receipt["rows"]
    samples = raw[1:-1:2]
    rbs = raw[2:-1:2]
    nsquared = n * n
    max_logZ = max(x["logZ"] for x in samples)
    unnormalized = []
    for rep, (x, y) in enumerate(zip(samples, rbs, strict=True)):
        assert x["kind"] == "replicate" and y["kind"] == "firstmerge_rb_calibration"
        assert x["replicate"] == y["replicate"] == rep and x["tau"] == y["tau"] == 2
        assert not x["dead"] and math.isfinite(x["logZ"])
        seed = (0xA022000000000000 + int(m["domain"], 16) + n * 1_000_000 + 20000 + rep) % (1 << 64)
        assert int(x["seed_hex"], 16) == seed
        hist = y["mean_m_hist"]
        assert len(hist) == 2 * n - 1 and hist[0] == hist[1] == 0
        assert all(math.isfinite(v) and v >= 0 for v in hist)
        assert y["mean_m"] >= 1 and math.isclose(sum(hist), y["mean_m"], rel_tol=1e-12)
        assert math.isclose(y["mean_m_F0"] - 2*y["mean_m_G0"], 4*y["mean_m_d2"], rel_tol=1e-12)
        unnormalized.append(math.exp(x["logZ"] - max_logZ))
    weights = [w*y["mean_m"] for w, y in zip(unnormalized, rbs)]
    Z = sum(weights)
    hist = [sum(w*y["mean_m_hist"][d] for w, y in zip(unnormalized, rbs))/Z
            for d in range(2*n - 1)]
    q_m = sum(w*y["mean_m_d2"] for w, y in zip(unnormalized, rbs))/(Z*nsquared)
    q_l = sum(w*x["tau2_lift_scaled_mean_L_distance2"] for w, x in zip(unnormalized, samples))
    q_l /= nsquared*sum(w*x["tau2_lift_scaled_mean_L"] for w, x in zip(unnormalized, samples))
    weight_n_eff = Z * Z / sum(w * w for w in weights)
    return dict(n=n, P=P, R=R, L=m["move_steps"], domain=m["domain"], Q_mRB=q_m, Q_L=q_l,
                far_n_over_2=sum(hist[n//2 + 1:]), distance_histogram_mRB=hist,
                normalized_rep_weight_max=max(weights)/Z, replicate_weight_ess=weight_n_eff,
                seed_hex=[x["seed_hex"] for x in samples], wall_seconds=raw[-1]["wall_seconds"],
                reported_gpu_seconds=sum(x["gpu_seconds"] for x in samples),
                diagnostic="finite SMC: no mixing or error bound is claimed")


if __name__ == "__main__":
    data = check(*(Path(x) for x in sys.argv[1:5]))
    print(json.dumps(data, indent=2))
