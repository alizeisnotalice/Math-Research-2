"""Generate six rational inputs and transfer the six-body area bound to them."""
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import sys

from flint import arb, ctx
from certify_large_support import curvature
from certify_nine import S, append_point, ar, arc_body, fourier_body

HERE = Path(__file__).resolve().parent
PI = F(22, 7)
DISK, ARC, FOURIER = 1408, 1024, 2048


def main(output):
    ctx.prec = 160
    pi = arb.pi()
    assert pi < ar(PI)
    data = json.loads((HERE / 'large-support-harmonic-input.json').read_bytes())
    terms = [(k, F(c)) for k, c in data['terms']]
    lower_radius = curvature(terms, pi)
    upper_radius = 1 - lower_radius
    # Width one implies rho(theta) + rho(theta+pi) = 1.
    assert all(k % 2 for k, c in terms)
    assert (F(1, 2) + sum(abs(c) for k, c in terms))**2 + sum(k*abs(c) for k, c in terms)**2 < 1
    error = max(PI**2/(4*DISK**2), PI**2/(8*ARC**2),
                upper_radius * PI**2/(2*FOURIER**2)) + F(12, S)
    target = F(json.loads((HERE / 'tomar-six-high-input.json').read_bytes())['target'])
    assert target == F('.837295006222')
    area = (target - PI*error**2)/(1 + 2*error/(F(1, 2)-error))
    assert area > F('.8372908')
    groups = []
    diameters = []
    for n in [0, 3, 5, 7, 9, 'S']:
        pts = []
        if n == 0:
            for j in range(DISK):
                sine, cosine = (2*pi*j/DISK).sin_cos()
                append_point(pts, cosine/2, sine/2)
        elif n == 'S':
            fourier_body(pts, FOURIER, pi, terms, arb(0), arb(0), arb(0))
        else:
            arc_body(pts, (ARC+n-1)//n, pi, n, arb(0), arb(0), arb(0))
        pts = sorted(set(pts))
        diameter = max((x-u)**2 + (y-v)**2 for (x,y),(u,v) in combinations(pts, 2))
        assert diameter <= (S+10)**2
        groups.append(dict(body=n, points=pts))
        diameters.append(diameter)
    raw = (json.dumps(dict(denominator=S+10, groups=groups), separators=(',', ':'))+'\n').encode()
    output.write_bytes(raw)
    print(json.dumps(dict(group_sizes=[len(g['points']) for g in groups], denominator=S+10,
                          diameter_squared_numerators=diameters, hausdorff_error=str(error),
                          area_lower=str(area), strict_floor='.8372908', sha256=sha256(raw).hexdigest()), indent=2))


if __name__ == '__main__':
    main(Path(sys.argv[1]))
