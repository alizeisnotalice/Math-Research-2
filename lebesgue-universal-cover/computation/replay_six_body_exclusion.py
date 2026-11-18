"""Replay the nine-harmonic exclusion with the original affine area verifier."""
from fractions import Fraction as F
from hashlib import sha256
import json
from math import lcm
from pathlib import Path
import sys
from flint import arb, ctx
import certify_five_harmonic as fan
from certify_large_support import curvature

def replay(path, witness):
    ctx.prec = 160
    pi = arb.pi()
    raw = path.read_bytes()
    spec = json.loads(raw)
    cert = json.loads(witness.read_bytes())
    assert sha256(raw).hexdigest() == cert["input_sha256"]
    terms = [(k, F(c)) for k, c in spec["terms"]]
    assert terms == [(k, F(c)) for k, c in cert["test_harmonics"]]
    assert curvature(terms, pi) == F(1, 100)
    fan.NEW_RADIUS = F(1, 2) - sum(abs(c) for _, c in terms)
    fan.T_LIMIT = F(18, 25) - fan.NEW_RADIUS
    fan.L = sum(k*abs(c) for k, c in terms)
    assert 0 < fan.NEW_RADIUS and 0 < fan.T_LIMIT
    assert pi < fan.ar(fan.PI_UPPER)
    for n, _, _, r, lip in fan.OLD:
        assert 1-1/(2*(pi/(2*n)).cos()) > fan.ar(r)
        assert (pi/(2*n)).tan()/2 < fan.ar(lip)
    center, half = (list(map(F, spec[k])) for k in ("center", "halfwidths"))
    assert len(center) == len(half) == 11 and min(half) > 0
    old, _ = fan.old_common_polygon(center, half, pi, k=cert["old_resolution"])
    material = []
    fan.fourier_body(material, cert["test_material_count"], pi, terms, arb(0), arb(0), arb(0))
    rows = cert["angular_witnesses"]
    den = lcm(*(F(endpoint).denominator for endpoint, _ in rows))
    converted = [(int(F(stop)*den), [[x, y, w] for (x,y), w in sites]) for stop, sites in rows]
    target = F(spec["target"])
    assert target == F(cert["target"])
    lower = fan.verify(old, material, pi, target, converted, den)
    radius = fan.ar(fan.NEW_RADIUS + fan.T_LIMIT)
    far = (pi-(1/(2*radius)).acos())/4+(radius**2-fan.ar(F(1,4))).sqrt()/2
    assert far > fan.ar(lower)
    return dict(input=path.name, cells=len(rows), lower=str(lower), target=str(target))

if __name__ == "__main__":
    print(json.dumps(replay(*(Path(s) for s in sys.argv[1:])), indent=2))
