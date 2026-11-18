"""Bound all positive-turn hinge cuts using the new compulsory fold endpoints."""
import argparse
from contextlib import redirect_stdout
from fractions import Fraction as F
from hashlib import sha256
import io
import json
from pathlib import Path
import runpy

HERE = Path(__file__).resolve().parent


def certificate():
    def read(name):
        return json.loads((HERE / name).read_bytes())

    def digest(name):
        return sha256((HERE / name).read_bytes()).hexdigest()

    fold_name = "ah-ideal-root-fold-certificate.json"
    walls = read("ah-compulsory-fold-walls-certificate.json")
    assert walls["fold_certificate_sha256"] == digest(fold_name)
    pivot = F("0.428684848443")
    lo, hi = map(F, read(fold_name)["unique_fold_exact_enclosure"])
    assert lo < hi < pivot < F(".429")

    # Replay the original Green formula and its earlier regression before
    # changing its translation origin to the newly forced straight endpoint.
    with redirect_stdout(io.StringIO()):
        old = runpy.run_path(str(HERE / "certify_ah_compulsory_chord_area.py"))
        fixed = runpy.run_path(str(HERE / "bound_fixed_hinge_family.py"))
    density = old["density"]
    g = density.__globals__
    I = g["I"]
    g["W"] = g["add"](g["z"], g["mul"](I(pivot), g["D"]))
    cells, order = old["cells"], old["order"]
    width = pivot / cells
    integral = I(0)
    for i in range(cells):
        mid = (i + F(1, 2)) * width
        coeff = density(I(mid))
        for k in range(0, order, 2):
            integral += width**(k+1) / ((k+1)*2**k) * coeff[k]
        remainder = density(I(i*width, (i+1)*width))[order]
        integral += width**(order+1) / ((order+1)*2**order) * remainder
    cap = integral / 2
    gain = fixed["gain_endpoint"](fixed["gain"].upper())
    floor = old["old_area"].lo - gain - cap.hi
    improvement = floor - fixed["floor"]
    assert cap.lo > 0 and improvement > 0
    return dict(pivot=str(pivot), fold_sha256=digest(fold_name),
                wall_sha256=digest("ah-compulsory-fold-walls-certificate.json"),
                cells=cells, order=order, cap=cap.out(), EH_gain_upper=str(gain),
                family_floor=str(floor), increase_over_previous=str(improvement),
                family_floor_down_15=f"0.{int(floor*10**15):015d}")


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("output", type=Path)
    args = cli.parse_args()
    args.output.write_text(json.dumps(certificate(), indent=2) + "\n")
