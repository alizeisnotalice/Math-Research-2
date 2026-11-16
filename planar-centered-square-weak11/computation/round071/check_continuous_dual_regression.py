#!/usr/bin/env python3
"""Actual deterministic exact regressions for the R71 continuous oracle."""
import argparse
import hashlib
import json
from pathlib import Path
import platform
import sys
import time
from continuous_dual import SCHEMA, check_dual


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--input", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    start = time.perf_counter()
    raw = a.input.read_bytes()
    cases = json.loads(raw)["cases"]
    results = []
    for fixture in cases:
        data = {k:v for k,v in fixture.items() if not k.startswith("expected_") and k != "name"}
        data["schema"] = SCHEMA
        result = check_dual(data)
        assert result["passed"] == fixture["expected_passed"], fixture["name"]
        for key in ("checked_cells", "minimum_slack"):
            if "expected_"+key in fixture:
                assert result[key] == fixture["expected_"+key], (fixture["name"], key)
        if "expected_minimum_dimension" in fixture:
            assert result["minimum_cell"]["dimension"] == fixture["expected_minimum_dimension"]
        # Independently evaluate every reported representative using a direct
        # membership path, avoiding the oracle's axis masks/signature cache.
        from fractions import Fraction as F
        for cell in result["cells"]:
            point = tuple(map(F,cell["representative"]))
            total = F(0)
            for sign, family in ((-1,data["witnesses"]),(1,data["caps"])):
                active = []
                for r in family:
                    x0,x1,y0,y1 = map(F,r["bounds"])
                    x,y = point
                    inside = (x0<=x<=x1 and y0<=y<=y1) if r["boundary"] == "closed" else (x0<x<x1 and y0<y<y1)
                    if inside:
                        active.append(r["id"])
                        total += sign*F(r["weight"])
                key = "active_witness_ids" if sign == -1 else "active_cap_ids"
                assert active == cell[key]
            assert total == F(cell["slack"])
        results.append({"name":fixture["name"],"passed_assertions":True,"result":result})
    # Fail closed for inexact numbers, negative weights and unauthorized
    # zero-capacity exemptions; these must never produce a valid certificate.
    bad = [
        {"schema":SCHEMA,"witnesses":[{"id":"w","bounds":[0,1,0,1],"weight":0.1,"boundary":"open"}]},
        {"schema":SCHEMA,"witnesses":[{"id":"w","bounds":[0,1,0,1],"weight":-1,"boundary":"open"}]},
        {"schema":SCHEMA,"zero_capacity_cap_ids":["absent"]},
        {"schema":SCHEMA,"caps":[{"id":"c","bounds":[0,1,0,1],"weight":1,"boundary":"open"}]}
    ]
    for data in bad:
        try:
            check_dual(data)
        except ValueError:
            pass
        else:
            raise AssertionError("invalid data accepted")
    out = {"schema":"r071-continuous-regression-result-v1","all_passed":True,
           "case_count":len(cases),"invalid_input_rejections":len(bad),"cases":results,
           "receipt":{"python":sys.version,"platform":platform.platform(),
                      "input_sha256":hashlib.sha256(raw).hexdigest(),
                      "script_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                      "oracle_sha256":hashlib.sha256(Path(__file__).with_name("continuous_dual.py").read_bytes()).hexdigest(),
                      "wall_seconds":time.perf_counter()-start}}
    a.output.write_text(json.dumps(out,indent=2)+"\n")
    print(json.dumps({"all_passed":True,"case_count":len(cases),"invalid_input_rejections":len(bad)}))


if __name__ == "__main__":
    main()
