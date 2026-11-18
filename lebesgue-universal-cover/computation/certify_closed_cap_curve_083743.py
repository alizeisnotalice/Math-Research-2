"""Bind the optimized cap comparison to the complete signed .83743 cover."""
import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path

from certify_cap_parameter_curve import verify as curve
from certify_refined_cap_window import verify as windows


def verify(comp, tree, production):
    def load(name):
        return json.loads((comp / name).read_bytes())

    terminal = load("finite-point-lift-083743.json")
    shrink = load("signed-shallow-cap-083743-replay.json")
    assert shrink["complete"] and shrink["terminal_sha256"] == sha256(
        (comp / "finite-point-lift-083743.json").read_bytes()).hexdigest()
    expected = {tuple(r[:2]): r for r in terminal["rows"]}
    observed = {}
    for name in ("lift83743", "remote265-collected", "remote1582-collected",
                 "remote688-collected", "remote16-collected", "remote-tail6-collected"):
        for path in (production / name).glob("part*.json"):
            if path.is_symlink():
                continue
            raw = json.loads(path.read_bytes())
            k = raw["part"], raw["pos"]
            st = raw["stat"]
            if not st[0]:
                continue
            row = [*k, raw["level"], *st[1:5], st[6], st[-1]]
            assert k not in observed or observed[k] == row
            observed[k] = row
    assert observed == expected
    row = curve(comp)["samples"][1]
    window = windows(comp, tree)
    g = F(row["integer_goal"])
    assert F(shrink["goal"]) == g == F(window["conditional_rows"][0]["integer_goal"])
    assert len(expected) == shrink["ordinary"] == 2572
    exceptions = {tuple(r[:2]) for r in shrink["high_rows"]}
    assert exceptions == set(map(tuple, terminal["five_only"])) and len(exceptions) == 87
    assert all(F(r[2])-F(window["exceptional_loss"]) > g for r in shrink["high_rows"])
    assert F(row["conditional_floor"]) > F(".83743127005895")
    return dict(ordinary=len(expected), ordinary_nodes=sum(r[3] for r in expected.values()),
                exceptional=len(exceptions), integer_goal=str(g), parameters=row,
                general_lower_display=".83743127005895")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ("tree", "production", "out"):
        parser.add_argument(key, type=Path)
    args = parser.parse_args()
    result = verify(Path(__file__).resolve().parent, args.tree, args.production)
    args.out.write_text(json.dumps(result, indent=2)+"\n")
