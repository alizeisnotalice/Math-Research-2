"""Explicit six-set catalogues for a host of area at most the selected ceiling."""
import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path

from certify_finite_083740_transfer import certified_core_scale
from certify_polygon_083731 import diameter_squared
from tangent_polygons_finite_cover import ar, arb, ctx, generate

COMP = Path(__file__).resolve().parent

def write_catalogue(ceiling, goal, out):
    groups, sampling = generate(ceiling, goal)
    core = min(certified_core_scale(g["body"], g["points"], g["denominator"])[0]
               for g in groups[:-1])
    assert 0 < core <= 1
    diameters = [diameter_squared(group) for group in groups]
    assert all(v <= 1 for v in diameters)
    area = core**2 * goal
    assert area > ceiling
    name = "obstacles-" + str(ceiling.numerator) + "-" + str(ceiling.denominator) + ".json"
    raw = (json.dumps(groups, separators=(",", ":")) + "\n").encode()
    (out / name).write_bytes(raw)
    # Compare with the necessary vertex count in the disk-containment transfer.
    m = sampling["counts"][0] - 1
    def best_possible(k):
        t = arb.pi()/k
        return ar(goal)*(t.cos()/(t/2).cos())**2
    while best_possible(m) >= ar(ceiling):
        m -= 1
    assert best_possible(m) < ar(ceiling) < best_possible(m+1)
    return dict(ceiling=str(ceiling), counts=sampling["counts"], core_scale=str(core),
                area_floor=str(area), squared_diameters=list(map(str, diameters)),
                necessary_disk_vertices=m+1, coordinate_file=name, sha256=sha256(raw).hexdigest())

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("out", type=Path)
    out = parser.parse_args().out
    ctx.prec = 160
    ledger = (COMP / "tomar-083742-stage-ledger.json").read_bytes()
    terminal = (COMP / "finite-point-lift-083742.json").read_bytes()
    data = json.loads(terminal)
    assert data["outer_audit_sha256"] == sha256(ledger).hexdigest()
    assert data["point_sha256"] == sha256((COMP / "finite-point-leaf-input.json").read_bytes()).hexdigest()
    goal = F(data["target"])
    assert goal == F(".83742") and F(json.loads(ledger)["integer_goal"], 2**57) > goal
    out.mkdir(parents=True, exist_ok=True)
    rows = [write_catalogue(F(t), goal, out) for t in (".83735", ".83740", ".83741")]
    print(json.dumps(dict(terminal_sha256=sha256(terminal).hexdigest(), goal=str(goal),
                          catalogues=rows), indent=2))

if __name__ == "__main__":
    main()
