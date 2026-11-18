"""Double and multiprecision numeric checks for the certified finite witness."""
import json
from pathlib import Path

from verify import Point, high_precision, integrate


def main():
    point = Point(1.25, 1.5, 0.125, 2**-8, 0.125, 1, 2**-40, 0, 4, 4)
    report = {"input": vars(point)}
    for name, baseline in (("corrected", False), ("matched_quartic", True)):
        report[name] = {
            "mpmath_65digits": high_precision(point, baseline),
            "double": integrate(point, baseline),
        }
    Path("witness-numerics.json").write_text(json.dumps(report, indent=2) + "\n")


if __name__ == "__main__":
    main()
