"""Directed bounds for the three regular constant-width inputs and the maximum-cut obstruction."""
from argparse import ArgumentParser
import json
from pathlib import Path
from flint import arb, ctx


def certify():
    ctx.prec = 192
    pi = arb.pi()
    t = arb(154)/10000
    threshold = arb(23373)/10000000
    jacobian = lambda n: 2*n*(pi/(2*n)).tan()-pi
    j3, j5, j7 = [jacobian(m) for m in (3, 5, 7)]
    cross = 2*(pi/30).tan()-pi/15
    d3, d5 = [(pi/(2*m)).tan()/2 for m in (3, 5)]
    l3, l5 = d3+t, d5+t
    l35 = l3+l5
    r5 = 1/(2*(pi/10).cos())
    assert t/l3 < pi/6 and t/l5 < pi/10
    assert (r5-arb(1)/2+2*t)/l35 < pi/6
    assert 3*t*t/l3 > threshold and 5*t*t/l5 > threshold

    # Fold the six zeros onto the 5-body's sign periods.
    count = 256
    a = pi/10
    sep = 2*a/3
    profile = lambda z: r5*(a-z).cos()-arb(1)/2
    energy_floor = arb(1)
    for i in range(count):
        low, high = a*i/(3*count), a*(i+1)/(3*count)
        vals = (profile(low), profile(sep+low), profile(sep-high))
        lower = 2*sum(max(arb(0), val.lower())**2 for val in vals)
        energy_floor = min(energy_floor, lower.lower())
    joint = energy_floor/(l35+l3+3*l5/5)
    assert joint > threshold
    base = pi/4+(2*j3+j5-cross)/16
    area = base+threshold/16
    assert area > arb("0.83253507411889576864")

    # Verify positivity and the random-cut domination for all m >= 7.
    a, b, d = (j5+cross)/2, (2*j3-j5-cross)/2, (j5-cross)/2
    cut_slack = 2-j7*(1/b+1/d)
    assert min(a, b, d, cut_slack) > 0
    assert (pi/4+(a+b+d)/8-base).contains(0)
    return {"regular_orders": [3, 5], "angular_intervals": count, "precision_bits": ctx.prec,
            "translation_threshold": "77/5000", "zero_energy_lower": str(energy_floor),
            "individual_penalty_3": str(3*t*t/l3), "individual_penalty_5": str(5*t*t/l5),
            "combined_penalty": str(joint), "rational_penalty": "23373/10000000",
            "maximum_cut_value": str(base), "cut_feasibility_margin": str(cut_slack),
            "area_lower": str(area), "strict_decimal_lower": "0.83253507411889576864"}


if __name__ == "__main__":
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    result = certify()
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(result["strict_decimal_lower"], result["combined_penalty"])
