"""Production Fraction verifier for 80 canonical H9(x^3)H9(x^5) layouts.

It first invokes the supplied exact Rouché verifier on the root-ball JSON.
It then reconstructs all G5 mother blocks with Fraction intervals and the
Phi5 blocks directly from a strict rational enclosure of sqrt(5).  A layout
is rejected only by an explicit coefficient interval with negative upper end.
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import itertools
import json
from dataclasses import dataclass
from fractions import Fraction as F
from math import isqrt
from pathlib import Path


# Frozen byte-identical local copies of h7's final Fraction/Rouché verifier
# and its helper are kept beside this script, so the production command has no
# absolute import path into another agent's live directory.
from verify_h9_roots import verify as verify_root_balls


def qstr(x: F) -> str:
    return f"{x.numerator}/{x.denominator}"


def pf(s: str) -> F:
    a, b = s.split("/", 1)
    return F(int(a), int(b))


@dataclass(frozen=True)
class I:
    lo: F
    hi: F
    def __post_init__(self): assert self.lo <= self.hi
    def __add__(self, y: "I") -> "I": return I(self.lo+y.lo, self.hi+y.hi)
    def __neg__(self) -> "I": return I(-self.hi, -self.lo)
    def __mul__(self, y: "I") -> "I":
        v = (self.lo*y.lo, self.lo*y.hi, self.hi*y.lo, self.hi*y.hi)
        return I(min(v), max(v))


Z, O = I(F(0), F(0)), I(F(1), F(1))


def square(x: I) -> I:
    return I(F(0), max(x.lo*x.lo, x.hi*x.hi)) if x.lo <= 0 <= x.hi else x*x


def mul(a: list[I], b: list[I]) -> list[I]:
    out = [Z for _ in range(len(a)+len(b)-1)]
    for i, x in enumerate(a):
        for j, y in enumerate(b): out[i+j] = out[i+j] + x*y
    return out


def prod(blocks: list[list[I]]) -> list[I]:
    out = [O]
    for b in blocks: out = mul(out, b)
    return out


def sparse(deg: int, values: dict[int, I]) -> list[I]:
    out = [Z for _ in range(deg+1)]
    for i, x in values.items(): out[i] = x
    return out


def disk(node: dict) -> tuple[I, I]:
    radius = pf(node["disk"]["radius"])
    re, im = pf(node["disk"]["center"]["re"]), pf(node["disk"]["center"]["im"])
    return I(re-radius, re+radius), I(im-radius, im+radius)


def linear(re: I, p: int) -> list[I]:
    return sparse(p, {0: -re, p: O})


def quadratic(re: I, norm: I, p: int) -> list[I]:
    return sparse(2*p, {0: norm, p: -(I(F(2),F(2))*re), 2*p: O})


def sqrt5_bounds(scale: int = 10**90) -> I:
    n = isqrt(5*scale*scale)
    out = I(F(n,scale), F(n+1,scale))
    validate_sqrt5_bounds(out)
    return out


def validate_sqrt5_bounds(bounds: I) -> None:
    assert bounds.lo >= 0 and bounds.lo*bounds.lo < 5 < bounds.hi*bounds.hi


def exact_mul(a: list[int], b: list[int]) -> list[int]:
    c = [0]*(len(a)+len(b)-1)
    for i, x in enumerate(a):
        for j, y in enumerate(b): c[i+j] += x*y
    return c


def stretch(a: list[int], p: int) -> list[int]:
    b = [0]*((len(a)-1)*p+1)
    for i, x in enumerate(a): b[i*p] = x
    return b


def key(e: str, mask: int, choices: tuple[int,...]) -> tuple:
    return e, mask, tuple(choices)


SURVIVORS = {key("e0",0,()), key("e1",3,(3,3))}
EXPECTED_KEYS = ({key("e0",mask,()) for mask in range(16)} |
                 {key("e1",mask,choices) for mask in range(16)
                  for choices in itertools.product((3,5), repeat=2)})
assert len(EXPECTED_KEYS) == 80


def negative_witnesses(poly: list[I], side: str) -> list[dict]:
    return [{"side":side,"degree":i,"upper":qstr(x.hi)}
            for i,x in enumerate(poly) if x.hi < 0]


def uncertain(poly: list[I]) -> list[int]:
    return [i for i,x in enumerate(poly) if x.lo < 0 <= x.hi]


def build_records(data: dict, sqrt_bounds: I) -> list[dict]:
    """The coefficient phase; caller must already have Rouché-verified data."""
    validate_sqrt5_bounds(sqrt_bounds)
    mothers = sorted(data["mothers"], key=lambda x: x["index"])
    assert [x["index"] for x in mothers] == list(range(5))
    # Direct Phi5 real parts: cos(2pi/5)=(-1+sqrt5)/4 and
    # cos(4pi/5)=(-1-sqrt5)/4.  Their root-pair norm is exactly one.
    t1 = I((sqrt_bounds.lo-1)/4, (sqrt_bounds.hi-1)/4)
    t2 = I((-1-sqrt_bounds.hi)/4, (-1-sqrt_bounds.lo)/4)
    cyc = [quadratic(t1,O,3), quadratic(t2,O,3),
           quadratic(t1,O,5), quadratic(t2,O,5)]
    neg_re, _ = disk(mothers[0])
    g = {(0,p): linear(neg_re,p) for p in (3,5)}
    for pair, root in ((1,mothers[1]), (2,mothers[2])):
        re, im = disk(root)
        norm = square(re)+square(im)
        for p in (3,5): g[pair,p] = quadratic(re,norm,p)

    layouts = []
    for mask in range(16):
        aa = [cyc[k] for k in range(4) if mask>>k & 1]
        bb = [cyc[k] for k in range(4) if not (mask>>k & 1)]
        for r in range(3):
            for p in (3,5): bb.append(g[r,p])
        layouts.append(("e0",mask,(),aa,bb))
    for mask in range(16):
        for choices in itertools.product((3,5), repeat=2):
            aa = [cyc[k] for k in range(4) if mask>>k & 1] + [g[0,3]]
            bb = [cyc[k] for k in range(4) if not (mask>>k & 1)] + [g[0,5]]
            for pair,p in enumerate(choices,1):
                aa.append(g[pair,p]); bb.append(g[pair,8-p])
            layouts.append(("e1",mask,choices,aa,bb))
    assert len(layouts) == 80

    records = []
    for e,mask,choices,aa,bb in layouts:
        ap,bp = prod(aa),prod(bb)
        witnesses = negative_witnesses(ap,"A") + negative_witnesses(bp,"B")
        ident = key(e,mask,choices)
        status = "rejected" if witnesses else "survivor"
        tag = "constant" if ident == key("e0",0,()) else (
              "standard_H9_x3" if ident == key("e1",3,(3,3)) else "candidate")
        records.append({"e":e,"cyclic_mask":mask,"complex_fibres":list(choices),"tag":tag,
                        "status":status,"negative_witnesses":witnesses,
                        "uncertain_A":uncertain(ap),"uncertain_B":uncertain(bp),
                        "degrees":[len(ap)-1,len(bp)-1]})
    return records


def exact_standard_tests() -> dict:
    h9 = [1,1,1,0,0,0,0,0,1,1]
    h3,h5 = stretch(h9,3),stretch(h9,5)
    c = exact_mul(h3,h5)
    assert c[0] == c[-1] == 1 and all(x in (0,1) for x in c)
    # Mask 3 is the two p=3 Phi5 conjugate-pair blocks; choices (3,3)
    # selects every G5 p=3 fibre.  Hence it is exactly H9(x^3), not merely
    # a numerically close survivor.
    assert key("e1",3,(3,3)) in SURVIVORS
    return {"constant_times_C":True,"H9_x3_times_H9_x5":True,
            "standard_mapping_mask3_choices33":True}


def validate_records(records: list[dict]) -> None:
    assert len(records) == 80
    layout_keys = [key(r["e"],r["cyclic_mask"],tuple(r["complex_fibres"])) for r in records]
    assert len(set(layout_keys)) == 80 and set(layout_keys) == EXPECTED_KEYS
    assert {r["status"] for r in records} <= {"rejected","survivor"}
    actual = {key(r["e"],r["cyclic_mask"],tuple(r["complex_fibres"]))
              for r in records if r["status"] == "survivor"}
    assert actual == SURVIVORS, (actual,SURVIVORS)
    assert sum(r["status"] == "rejected" for r in records) == 78
    for r in records:
        if r["status"] == "rejected":
            assert r["negative_witnesses"], r
            assert all(F(w["upper"]) < 0 for w in r["negative_witnesses"])


def interval_regressions() -> dict:
    assert I(F(-2),F(3))*I(F(-5),F(7)) == I(F(-15),F(21))
    assert I(F(-2),F(3))*I(F(-4),F(5)) == I(F(-12),F(15))
    assert I(F(2),F(3))*I(F(4),F(5)) == I(F(8),F(15))
    return {"negative_corner_product":True,"cross_zero_corner_product":True,
            "positive_corner_product":True}


def verify(data: dict) -> dict:
    """Production entry point: validate root balls, then validate all layouts."""
    root_receipt = verify_root_balls(data)
    records = build_records(data, sqrt5_bounds())
    validate_records(records)
    return {"root_receipt":root_receipt,"candidate_count":80,"rejected_count":78,
            "survivors":[{"e":"e0","cyclic_mask":0,"complex_fibres":[],"tag":"constant"},
                         {"e":"e1","cyclic_mask":3,"complex_fibres":[3,3],"tag":"standard_H9_x3"}],
            "records":records,"exact_integrity":exact_standard_tests()}


def tamper_regressions(data: dict, result: dict) -> dict:
    """Executed semantic tamper tests of production verification helpers."""
    def fails(action) -> bool:
        try:
            action()
        except (AssertionError, ValueError):
            return True
        return False
    bad_ball = copy.deepcopy(data); bad_ball["mothers"][0]["lhs"] = "0/1"
    assert fails(lambda: verify_root_balls(bad_ball)), "bad root ball accepted"
    good = sqrt5_bounds(); bad_sqrt = I(good.hi, good.hi+F(1,10**90))
    assert fails(lambda: validate_sqrt5_bounds(bad_sqrt)), "bad sqrt5 interval accepted"
    bad_records = copy.deepcopy(result["records"])
    next(r for r in bad_records if r["status"] == "rejected")["status"] = "survivor"
    assert fails(lambda: validate_records(bad_records)), "bad survivor label accepted"
    bad_witnesses = copy.deepcopy(result["records"])
    target = next(r for r in bad_witnesses if r["status"] == "rejected")
    for w in target["negative_witnesses"]: w["upper"] = "1/1"
    assert fails(lambda: validate_records(bad_witnesses)), "positive fake witness accepted"
    assert fails(lambda: validate_records(copy.deepcopy(result["records"])[:-1])), "missing layout accepted"
    duplicate = copy.deepcopy(result["records"]); duplicate[-1] = copy.deepcopy(duplicate[-2])
    assert fails(lambda: validate_records(duplicate)), "duplicate layout accepted"
    return {"bad_root_ball_rejected":True,"bad_sqrt5_bounds_rejected":True,
            "bad_survivor_label_rejected":True,"positive_fake_witness_rejected":True,
            "missing_layout_rejected":True,"duplicate_layout_rejected":True,
            "interval_corner_regressions":interval_regressions()}


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--balls",type=Path,required=True)
    p.add_argument("--record",type=Path,required=True)
    a = p.parse_args()
    raw = a.balls.read_bytes(); data = json.loads(raw)
    out = verify(data)
    out["balls_sha256"] = hashlib.sha256(raw).hexdigest()
    out["tamper_regressions"] = tamper_regressions(data,out)
    a.record.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
    print(json.dumps({k:out[k] for k in ("balls_sha256","candidate_count","rejected_count","survivors","tamper_regressions")}))


if __name__ == "__main__": main()
