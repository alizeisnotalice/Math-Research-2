"""Certify necessary radial translation bounds for D,B3,B5,B7,B9."""

from pathlib import Path
from hashlib import sha256
from fractions import Fraction as F
import json
import numpy as np
from scipy.spatial import ConvexHull
from flint import arb, ctx

ctx.prec = 160
w = Path(__file__).resolve().parent
f = w / "five-body-radial-material.npz"
tab = np.load(f)
Q = 2**28
pts = tab["points"]
ct = tab["counts"]
assert (
    pts.dtype == np.int32
    and ct.tolist() == [96, 96, 160, 224, 288]
    and pts.shape == (sum(ct), 2)
)
pts = pts.astype(np.int64)
ofs = np.r_[0, np.cumsum(ct)]


def a(q):
    v = F(q)
    return arb(v.numerator) / v.denominator


R = [".06835", ".0665", ".06004", ".05699"]
N = 2048


def area2(p, fixed, error):
    indices = ConvexHull(p.astype(float)).vertices
    v = p[indices]
    moving = (indices >= fixed).astype(np.int64)
    u = np.roll(v, -1, axis=0)
    next_moving = np.roll(moving, -1)
    dets = v[:,0]*u[:,1] - v[:,1]*u[:,0]
    # det(v+m*d,u+n*d) is affine in the common translation correction d.
    cx = moving*u[:,1] - next_moving*v[:,1]
    cy = next_moving*v[:,0] - moving*u[:,0]
    assert np.max(np.abs(v)) < Q and error < Q//10 and len(v) < 1500
    assert np.min(dets-error*(abs(cx)+abs(cy))) > 0
    assert np.count_nonzero((v[:,1]<=0)&(u[:,1]>0)) == 1
    return sum(map(int,dets))-error*(abs(int(sum(cx)))+abs(int(sum(cy))))

disc = pts[: ofs[1]]
assert all(int(x) ** 2 + int(y) ** 2 <= (Q // 2) ** 2 for x, y in disc)
out = []
for k, (n, rr) in enumerate(zip((3, 5, 7, 9), R), 1):
    rad = a(rr)
    vis = pts[ofs[k] : ofs[k + 1]]
    Rc = 1 / (2 * (arb.pi() / (2 * n)).cos())
    vertices = []
    for j in range(n):
        s, c = (2 * arb.pi() * j / n).sin_cos()
        vertices.append((Rc * c, Rc * s))
    for x, y in vis:
        ux, uy = arb(int(x)) / Q, arb(int(y)) / Q
        assert all((ux - X) ** 2 + (uy - Y) ** 2 <= 1 for X, Y in vertices)
    delta = rad*arb.pi()/(2*n*N)+arb(2)/Q
    error = int(float(delta*Q))+2
    assert delta*Q < error
    m = None
    for j in range(N):
        t = arb.pi()*(2*j+1)/(2*n*N)
        s,c = t.sin_cos()
        xr,yr = rad*c*Q,rad*s*Q
        tx,ty = round(float(xr.mid())),round(float(yr.mid()))
        assert abs(xr-tx)<1 and abs(yr-ty)<1
        twice=area2(np.r_[disc,vis+[tx,ty]],len(disc),error)
        if m is None or twice<m:m=twice
    lower=arb(m)/(2*Q**2)

    assert lower > a(".837295006222")
    out.append(
        dict(
            n=n,
            radius=rr,
            uniform_twice_area_numerator=m,
            coordinate_error_grid_units=error,
            area_lower=str(lower),
        )
    )
print(
    json.dumps(
        dict(
            source_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),
            source_table_sha256="fe5b816060d713642dc44e3f8981f4ae5d1b2b1c1e271b7dff75f031dd75d87c",
            target="0.837295006222",
            grid=N,
            points_sha=sha256(pts.tobytes()).hexdigest(),
            table_sha=sha256(f.read_bytes()).hexdigest(),
            rows=out,
        ),
        indent=2,
    )
)
