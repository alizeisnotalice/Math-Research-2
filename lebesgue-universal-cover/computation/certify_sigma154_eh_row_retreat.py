"""Recheck AH placements with a specified retreat of only the EH support row."""
from argparse import ArgumentParser
from fractions import Fraction as F
from hashlib import sha256
import inspect
import json
from pathlib import Path
import sys


def run(comp, inp, out):
    sys.path.insert(0, str(comp))
    from certify_sigma_extended_graph import geometry as ordinary

    cfg = json.loads(inp.read_bytes())
    retreat = F(cfg["EH_margin_retreat"])
    assert F(cfg["sigma_degrees"]) == F("1.54") and retreat > 0
    if "slope_bins" in cfg:
        from certify_sigma_derivative_brackets import run as solver
    else:
        from certify_sigma_jet_bracket import run as solver
    checked = {}

    def geometry(*a, **kw):
        ns, source, tail = ordinary(*a, **kw)
        assert F(ns["eh_depth"]) == F(".000555")
        A = ns["A"]
        rows = ns["rows"]
        selected, = [i for i,(label,_,_) in enumerate(rows) if label.startswith("U1-Te")]
        label, normal, offset = rows[selected]
        rows[selected] = label, normal, offset-A(retreat)
        checked["EH_tightening"] = dict(label=label, retreat=str(retreat), depth=ns["eh_depth"])
        checked["geometry_sha256"] = sha256(source.encode()).hexdigest()
        return ns, source, tail

    env = {**solver.__globals__, "geometry": geometry}
    exec(compile(inspect.getsource(solver), "retreat_one_EH_row", "exec"), env)
    env["run"](comp, inp, out)
    record = json.loads(out.read_bytes())
    assert record["parameters"] == cfg
    assert record.get("geometry_sha256",checked["geometry_sha256"]) == checked["geometry_sha256"]
    record.update(checked)
    record["source_hashes"][Path(__file__).name] = sha256(Path(__file__).read_bytes()).hexdigest()
    out.write_text(json.dumps(record,indent=2)+"\n")


if __name__ == "__main__":
    assert __debug__
    ap = ArgumentParser(description=__doc__)
    for k in ("comp","input","output"):
        ap.add_argument(k,type=Path)
    a = ap.parse_args()
    run(a.comp.resolve(),a.input,a.output)
