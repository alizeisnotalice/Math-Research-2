"""Certify area bounds and countermodels for the Round 59/60 width relaxations."""
import json
from fractions import Fraction as F
import platform
from flint import arb, ctx

N = 256
A = F(1, 32)
SHIFTS = [F(k, 1000) for k in (1,2,3,4,5,6,7,8,10,12)]


def ar(value):
    z = F(value)
    return arb(z.numerator)/z.denominator


def perturbed_weights(pi):
    rows = []
    for shift in SHIFTS:
        angles = [ar(shift), pi/4, pi/2, 3*pi/4]
        weights = []
        for j, theta in enumerate(angles):
            d = arb(1)
            for k, phi in enumerate(angles):
                if k != j:
                    d *= abs((theta-phi).sin())
            assert d > 0
            weights.append(1/d)
        rows.append((angles, weights, sum(weights)))
    return rows


def imposed_sup_upper(mu0, pi, rows):
    # Assume the excess-width mean <= A+mu0.  For each frequency use
    # the best among the exact equal-angle multiplier and ten perturbed ones.
    bounds = []
    for n in range(1,N+1):
        candidates = [ar(mu0).upper()] if n % 4 == 0 else []
        for angles, weights, total in rows:
            multiplier_x = sum(v*(2*n*a).cos() for v,a in zip(weights,angles,strict=True))
            multiplier_y = sum(v*(2*n*a).sin() for v,a in zip(weights,angles,strict=True))
            multiplier = (multiplier_x**2+multiplier_y**2).sqrt()
            if multiplier > 0:
                numerator = ar(A+mu0)*total-ar(F(1,4))
                assert numerator>0
                candidates.append((numerator/multiplier).upper())
        bounds.append(min(candidates))
    return ar(A+mu0)+2*sum(bounds)+ar((1+A+mu0)/(2*N+1))


def main():
    ctx.prec=160
    pi=arb.pi()
    equal59=(pi/8).sin()/8
    equal60=(arb(2).sqrt()+2*(pi/12).cos())*(pi/24).cos()/(arb(3).sqrt()*(2+arb(2).sqrt()))-arb(1)/2
    checks=[]
    rows=perturbed_weights(pi)
    for delta,condition,corner in ((F(1,400000),"r59",equal59),(F(1,125000),"r60",equal60)):
        sup=imposed_sup_upper(delta,pi,rows)
        assert sup<corner
        checks.append({"angular_max_condition":condition,
                       "forced_mean_excess_over_1_over_32":str(delta),
                       "contradicted_upper_on_max_width_excess":str(sup),
                       "area_lower":str((1+ar(A+delta))*pi/4)})
    lam=F(1001,1000)
    wmin=ar(lam)*(1+1/(pi/9).cos())/2
    w4=ar(lam)*(pi/18).cos()*(pi/72).cos()/(pi/9).cos()
    circumscribed_regular_nonagon=ar(lam**2*F(9,4))*(pi/9).tan()
    assert wmin>ar(F(33,32)) and w4>1+equal59
    smooth_centered_hull=31*pi/128+arb(1)/16
    triangle_centered_hull=pi/8+arb(3).sqrt()/4
    assert circumscribed_regular_nonagon<smooth_centered_hull<triangle_centered_hull<ar(F("0.8361"))
    print(json.dumps({
        "python":platform.python_version(),
        "arb_precision":ctx.prec,
        "fourier_frequencies":N,
        "radian_shifts":[str(s) for s in SHIFTS],
        "mean_only_area_lower":str(33*pi/128),
        "strict_improvements":checks,
        "scaled_regular_nonagon_model":{"scale":str(lam),"area":str(circumscribed_regular_nonagon),
            "minimum_width":str(wmin),"minimax_equally_spaced_four_widths":str(w4)},
        "coherent_round59_exact_minimum":str(smooth_centered_hull),
        "combined_coherent_round59_round60_exact_minimum":str(triangle_centered_hull),
    },indent=2))


if __name__ == "__main__":
    main()
