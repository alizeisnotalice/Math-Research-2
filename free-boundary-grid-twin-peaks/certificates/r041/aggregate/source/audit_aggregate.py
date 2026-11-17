#!/usr/bin/env python3
"""Exact, standard-library audit for r071 aggregate outputs and n=5 targets."""
import argparse, json
from fractions import Fraction
P=argparse.ArgumentParser();P.add_argument("output",nargs="?");P.add_argument("--r069");P.add_argument("--n5-moment");P.add_argument("--n5-merge");a=P.parse_args()
def load(p):
    with open(p) as f:return json.load(f)
result={}
if a.n5_moment or a.n5_merge:
    if not(a.n5_moment and a.n5_merge):raise SystemExit("both n5 references required")
    m,g=load(a.n5_moment),load(a.n5_merge)
    expected={"A1":88366736882654697600,"A2":7063261481362958280000,"D2":67631454448823894520000}
    got={"A1":int(m["A1"]),"A2":int(m["A2"]),"D2":int(m["squared_euclidean_sum"])}
    merge={"A2":int(g["A2_from_W"]),"D2":int(g["D_euclidean_squared_from_W"])}
    assert got==expected and merge=={k:expected[k] for k in merge}
    result["n5_expected_owner_crosscheck"]={"passed":True,"expected":expected,"Q":str(Fraction(expected["D2"],25*expected["A2"]))}
if a.output:
    d=load(a.output);assert d["workers"]==1 and d["manhattan_histogram_included"] is False
    rows=d["ordered_q_lowb"];A2=int(d["A2"]);D2=int(d["D_euclidean_squared_sum"])
    assert sum(int(x["count"]) for x in rows)==A2 and sum(int(x["D_euclidean_squared_sum"]) for x in rows)==D2
    result["output"]={"n":d["n"],"A1":int(d["A1"]),"A2":A2,"D2":D2,"Q":str(Fraction(D2,d["n"]*d["n"]*A2)),"table_exact":True}
    if a.r069:
        ref=load(a.r069);assert d["n"]==ref["n"]
        def tab(x):return {(r["q"],r["low_b"]):(int(r["count"]),int(r["D_euclidean_squared_sum"])) for r in x["ordered_q_lowb"]}
        assert (d["A1"],d["A2"],d["D_euclidean_squared_sum"],tab(d))==(ref["A1"],ref["A2"],ref["D_euclidean_squared_sum"],tab(ref))
        result["r069_aggregate_and_each_qb_passed"]=True
print(json.dumps(result,sort_keys=True,indent=2))
