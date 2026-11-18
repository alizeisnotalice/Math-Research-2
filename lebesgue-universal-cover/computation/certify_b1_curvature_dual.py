"""Continuous support dual for a deep cut and either reflected patch tangent."""
import argparse
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path

from flint import arb, ctx
from certify_ah_released_tangent import load_geometry

HERE = Path(__file__).resolve().parent


def positive_area(a, b, left, right):
    """Integrate the positive part of a sin(theta)+b cos(theta) on a short arc."""
    assert 0 < right-left < arb.pi()
    def values(angle):
        s, c = angle.sin_cos()
        return a*s+b*c, -a*c+b*s

    kl, pl = values(left)
    kr, pr = values(right)
    radius = (a*a+b*b).sqrt()
    candidates = []
    # At most one interior zero. A downward crossing gives a maximum of
    # the primitive, an upward crossing gives a minimum.
    if not (kl < 0 or kr < 0):
        candidates.append(pr-pl)
    if not (kl > 0 or kr > 0):
        candidates.append(arb(0))
    if not (kl < 0 or kr > 0):
        candidates.append(radius-pl)
    if not (kl > 0 or kr < 0):
        candidates.append(pr+radius)
    return max(v.upper() for v in candidates)


def verify():
    geometry = load_geometry(HERE)
    ctx.prec = 192
    A = geometry["A"]
    pi = arb.pi()
    phi = 61*pi/1800
    dot, reflect = geometry["dot"], geometry["reflection"]
    cut_angle = A(Fraction(".0175"))
    cut = geometry["mul"](-1, geometry["perp"](geometry["rotate"](
        geometry["direct"], (cut_angle.cos(), cut_angle.sin()))))
    pivot = geometry["add"](geometry["z"], geometry["mul"](A(Fraction(".430")), geometry["D"]))
    input_path = HERE / "b1-curvature-dual-weights.json"
    results = []
    def unit(x):
        return x.cos(), x.sin()
    for row in json.loads(input_path.read_bytes()):
        weights = [A(Fraction(v)) for v in row["weights"]]
        wa, wb, wf = weights
        assert min(weights) > 0
        nu = [A(Fraction(v)) for v in row["equalities"]]
        p, dp = geometry["point_jets"](A(Fraction(row["tangent"])))
        q, dq = reflect(p[2]), reflect(dp[2])
        tangent = geometry["mul"](-1/geometry["norm"](dq), geometry["perp"](dq))
        v = wf/(3*phi).sin()
        c = -wa*dot(cut,pivot)-wb*dot(tangent,q)+v/2+(nu[0]+nu[1])/2+nu[2]
        slope = arb(0)
        terms = [(cut,wa),(reflect(tangent),wb),(unit(phi),-v),
                 (unit(phi+pi/3),v),(unit(phi+2*pi/3),-v),
                 (unit(pi/3),-nu[0]),(unit(2*pi/3),-nu[1])]
        events = []
        for u, weight in terms:
            theta = arb.atan2(u[1],u[0])
            if theta < 0:
                theta += 2*pi
            sign = 1
            if theta > pi:
                sign = -1
                theta -= pi
            assert 0 < theta < pi
            c += weight*(1-sign+u[0])/2
            slope += weight*u[1]
            events.append((theta,sign*weight))
        events.sort(key=lambda term: float(term[0].mid()))
        ends = [arb(0),*(t for t,_ in events),pi]
        integral = arb(0)
        for i,(left,right) in enumerate(zip(ends,ends[1:])):
            a, b = -nu[2], -nu[3]
            for theta, weight in events[i:]:
                a -= weight*theta.cos()
                b += weight*theta.sin()
            integral += positive_area(a,b,left,right)
        bound = (c+abs(slope)+integral).upper()
        assert bound < 0
        results.append(dict(tangent=row["tangent"], upper=str(bound),
                            positive_integral_upper=str(integral.upper())))
    return dict(rows=results, precision=ctx.prec, input_sha256=sha256(input_path.read_bytes()).hexdigest())


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("out", type=Path)
    args = parser.parse_args()
    args.out.write_text(json.dumps(verify(), indent=2)+"\n")
