"""Integer grid sizes sufficient for the two contact-based universal area bounds."""
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
import sys

from flint import arb, ctx


def exact(q):
    return arb(q.numerator)/q.denominator


def support_error(n):
    c = (arb.pi()/(2*n)).cos()
    return 2*(1-c)/((3+c*c).sqrt()+1+c)


def run(source):
    ctx.prec = 192
    receipt = json.loads(source.read_bytes())
    area = exact(Fraction(receipt["host_area_decimal_upper"]))
    budget = exact(Fraction(receipt["B"]))
    rows = []
    for divisor in (1, 2):
        def cost(n):
            radius = support_error(n)/divisor
            return (1+4*radius)*area + arb.pi()*radius**2
        low, high = 2, 10000
        assert cost(low) > budget and cost(high) < budget
        while high-low > 1:
            mid = (high+low)//2
            if cost(mid) < budget:
                high = mid
            else:
                assert cost(mid) > budget
                low = mid
        rows.append({"error_divisor": divisor, "first_n": high,
                     "radius": str(support_error(high)/divisor),
                     "excess_at_previous": str(cost(low)-budget),
                     "remaining_at_first": str(budget-cost(high))})
    return {"bits": ctx.prec, "area": receipt["host_area_decimal_upper"], "B": receipt["B"],
            "area_input_sha256": sha256(source.read_bytes()).hexdigest(), "rows": rows}


if __name__ == "__main__":
    assert __debug__
    result = run(Path(sys.argv[1]))
    Path(sys.argv[2]).write_text(json.dumps(result, indent=2)+"\n")
    print([(r["error_divisor"], r["first_n"]) for r in result["rows"]])
