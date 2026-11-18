"""Evaluate the fixed paired-arc AH area obstruction at the 1.54-degree Gibbs angle."""

import argparse
from fractions import Fraction
from pathlib import Path
import sys


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    if not __debug__:
        raise RuntimeError("Use Python with assertions enabled.")
    computation = Path(__file__).resolve().parent
    sys.path.insert(0, str(computation))
    import certify_low_sigma_ah_budget as budget

    budget.CASES["1.54"] = {
        "force": ".429523",
        "fold": (".4295222514414", ".4295222514415"),
    }
    budget.TURN = Fraction(".018675")
    budget.run(computation, "1.54", args.output)


if __name__ == "__main__":
    main()
