"""Independently validate selected candidates with adaptive Gauss–Kronrod quadrature."""
import argparse
import json
import math
import warnings
from pathlib import Path

from scipy.integrate import IntegrationWarning, quad


def check(case, maximum=False, radius=.25):
    t, m, eps, delta = (case[key] for key in ("t", "m", "eps", "delta"))
    c = math.sqrt((5913*t**3+954*t**2+96*t-8)/(21*t+6))
    d, x = t/c, (3*t-2)/c
    log_target = -.5*(math.log(2*math.pi*(m+eps))+(x+delta)**2/(m+eps))
    means, errors, splits = [], [], []
    for side, weight, power in [(-1., c, 1), (1., d, 3)]:
        cap = min(radius if not maximum else math.inf,
                  math.nextafter((m/weight)**(1/power), 0.)*(1.-2e-15))
        assert weight*cap**power <= m
        center = max(0., side*delta)
        width = math.sqrt(eps+weight*center**power)
        scales = [math.sqrt(eps), (eps/weight)**(1/power),
                  weight if power == 1 else 1/weight]
        knots = {0., cap, center}
        for i in range(-12, 13):
            fac = 2.**i
            knots.update([center+width*fac, center-width*fac])
            knots.update(s*fac for s in scales)
        for i in range(-12, 13):
            knots.add(-side*x+i*.5)
        knots = sorted({math.sqrt(k) for k in knots if 0 <= k <= cap})
        def normalized_integrand(s):
            u = s*s
            v = eps+weight*u**power
            cost = ((x+side*u)**2+(delta-side*u)**2/v)*.5
            return 2*s*math.exp(-cost-.5*math.log(v)-math.log(2*math.pi)-log_target)
        with warnings.catch_warnings():
            warnings.simplefilter("error", IntegrationWarning)
            estimates = [quad(normalized_integrand, a, b, epsabs=1e-15, epsrel=2e-12,
                              limit=100) for a,b in zip(knots[:-1],knots[1:],strict=True)]
        means.append(math.fsum(v for v,e in estimates))
        errors.append(math.fsum(e for v,e in estimates))
        splits.append(len(estimates))
    return dict(case=case, contributions=means, error_estimates=errors, split_counts=splits,
                ratio=math.fsum(means), independent_log_ratio=math.log(math.fsum(means)))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("runs", nargs="+", type=Path)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    rows = []
    for run in args.runs:
        reports = [json.loads(p.read_text()) for p in run.glob("shard-*.json")]
        assert len(reports) == 15
        for mode in ("radius_capped", "variance_capped"):
            cases = sorted(
                [item for r in reports for item in r["numerical_summary"][mode]["top"]],
                key=lambda item: item["log_ratio_order_96"], reverse=True,
            )[:5]
            rows.extend(dict(run=run.name, mode=mode, **check(case, mode == "variance_capped"))
                        for case in cases)
    for m in [1.125, 1.25, 1.5, 2., 12.]:
        rows.append(dict(mode="simple_parameters",
                         **check(dict(t=1./3, m=m, eps=2.**-20, delta=1./8))))
    args.out.write_text(json.dumps(rows, indent=2)+"\n")
    for row in rows:
        print(row.get("run"), row["mode"], row["case"]["m"], row["ratio"])


if __name__ == "__main__":
    main()
