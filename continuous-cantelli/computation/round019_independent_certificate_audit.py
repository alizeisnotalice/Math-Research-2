"""Independent, non-mutating arithmetic audit of the r019 q=1000 witness."""
from fractions import Fraction as Q
from math import factorial
from pathlib import Path
import hashlib
import json
import platform

import mpmath as mp


SOURCE_ROOT = Path(__file__).resolve().parent
SOURCE = SOURCE_ROOT / "round019_critical_tail_certificate.py"
SOURCE_RESULT = SOURCE_ROOT / "round019-critical-results.json"
OWN_RESULT = SOURCE_ROOT / "round019-independent-certificate-results.json"


def main() -> None:
    q = 1000
    P = Q(1, 2 * 3**9) * Q(4, 11) ** 1415
    B = Q(3, 8) ** 1500

    # Strict elementary inputs used by the paper proof.
    # e<11/4: sum through degree 2 is 5/2 and the factorial tail is <1/4.
    assert Q(5, 2) + Q(1, 4) == Q(11, 4)
    # e>8/3 from the degree-3 partial sum; cos(1)>1/2 from an alternating
    # lower partial sum; pi<22/7 entails sqrt(2*pi)<3.
    assert Q(1) + Q(1) + Q(1, 2) + Q(1, 6) == Q(8, 3)
    assert Q(1) - Q(1, 2) + Q(1, 24) - Q(1, 720) > Q(1, 2)
    assert Q(44, 7) < 9
    # e^7>1000, hence log(1000)<7; e^700>2^700>1000 fixes the phase.
    assert sum(Q(7**j, factorial(j)) for j in range(16)) > 1000
    assert 2**700 > 1000
    assert Q(q * q, 14) > 1500 and Q(q * q, 8) > 1500
    assert Q(9900, 7) < 1415

    # Exact endgame.  Positive annulus P, negative budget 4B, target B.
    assert P > 8 * B
    assert P - 4 * B > P / 2 > 4 * B > B

    source_json = json.loads(SOURCE_RESULT.read_text())
    assert source_json["exact_comparisons"]["positive_lower"] == str(P)
    assert source_json["exact_comparisons"]["total_negative_upper"] == str(4 * B)
    assert source_json["exact_comparisons"]["gaussian_modulus_upper"] == str(B)
    mp.mp.dps = 60
    ratio = mp.mpf(P.numerator) / P.denominator / (mp.mpf(B.numerator) / B.denominator)
    result = {
        "scope": "independent exact arithmetic audit; analytic interval bounds documented in report",
        "q": q,
        "P_over_B": mp.nstr(ratio, 50),
        "P_gt_8B": True,
        "rotated_lower_gt_B": True,
        "source_script_sha256": hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
        "source_json_sha256": hashlib.sha256(SOURCE_RESULT.read_bytes()).hexdigest(),
        "python": platform.python_version(),
        "mpmath": mp.__version__,
    }
    OWN_RESULT.write_text(json.dumps(result, indent=2) + "\n")
    print("r019 critical certificate independent arithmetic audit passed")
    print("P/B =", result["P_over_B"])


if __name__ == "__main__":
    main()
