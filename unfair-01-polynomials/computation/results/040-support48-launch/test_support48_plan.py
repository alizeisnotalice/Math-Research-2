#!/usr/bin/env python3
"""Local independent enumeration/count audit plus two real chunk computations."""
from __future__ import annotations
import copy,gzip,hashlib,itertools,json,math,sys,time
from pathlib import Path

HERE=Path(__file__).resolve().parent
def independent_supports(d,cap=6):
    for size in range(2,min(d+1,cap)+1):
        for inside in itertools.combinations(range(1,d),size-2):yield (0,*inside,d)
def count(d):return sum(math.comb(d-1,k) for k in range(min(4,d-1)+1))
def reject(label,f):
    try:f()
    except (RuntimeError,AssertionError,KeyError,IndexError):return label
    raise AssertionError("accepted mutation "+label)
def main():
    began=time.monotonic(); plan=json.loads((HERE/"plan/support48-plan.json").read_text())
    # Do not import support48_prepare: independently reconstruct the 47 pairs.
    expected=[];coverage=[]
    for total in (47,48):
        for m in range(1,total//2+1):
            ss=list(independent_supports(m));coverage.append({"total":total,"m":m,"n":total-m,"expected_I_count":len(ss),"binomial_count":count(m)})
            expected += [{"total":total,"m":m,"n":total-m,"index":i,"I":list(I)} for i,I in enumerate(ss)]
    assert coverage==plan["coverage"] and expected==plan["jobs"] and len(expected)==100005
    # Execute two concrete chunks locally using the exact frozen source and
    # remote worker body.  This is a local smoke, never SSH/production output.
    frozen=(HERE/"plan/frozen_support.py").read_text().split("\ndef main():",1)[0]
    env={"SOURCE_SHA":hashlib.sha256((HERE/"plan/frozen_support.py").read_bytes()).hexdigest()}
    exec(frozen+"\n"+(HERE/"remote_support48.py").read_text(),env)
    jobs=[next(x for x in plan["jobs"] if (x["total"],x["m"],x["index"])==(47,23,0)),next(x for x in plan["jobs"] if (x["total"],x["m"],x["index"])==(48,24,10902))]
    smoke=[]
    for job in jobs:
        row=env["compute_support48"](job); packed=__import__("base64").b64decode(row["gzip_base64"]);raw=gzip.decompress(packed)
        assert hashlib.sha256(packed).hexdigest()==row["gzip_sha256"] and hashlib.sha256(raw).hexdigest()==row["raw_sha256"]
        assert row["counts"]["raw"]==count(job["n"])
        smoke.append({"key":[job["total"],job["m"],job["n"],job["index"]],"raw_supports":row["counts"]["raw"],"records":len(raw.splitlines()),"gzip_sha256":row["gzip_sha256"]})
    # Use the package validator for its second, full job-set validation only
    # after the independent coverage equality above has passed.
    from support48_prepare import validate
    duplicate=copy.deepcopy(plan);duplicate["jobs"].append(copy.deepcopy(duplicate["jobs"][0]))
    wrong_i=copy.deepcopy(plan);wrong_i["jobs"][0]["I"]=[0,99]
    out_of_range=copy.deepcopy(plan);out_of_range["jobs"][0]["total"]=46
    rejected=[reject("duplicate_job",lambda:validate(duplicate)),reject("wrong_I",lambda:validate(wrong_i)),reject("out_of_range_total",lambda:validate(out_of_range))]
    result={"status":"passed_local_two_job_smoke_and_independent_counts","scope":"Local only; no SSH/MacBook production invocation.","jobs":len(expected),"degree_pairs":len(coverage),"independent_binomial_total":sum(x["binomial_count"] for x in coverage),"smoke":smoke,"negative_tests_rejected":rejected,"elapsed_seconds":time.monotonic()-began,"sha256":{"plan":hashlib.sha256((HERE/"plan/support48-plan.json").read_bytes()).hexdigest(),"frozen":env["SOURCE_SHA"],"remote":hashlib.sha256((HERE/"remote_support48.py").read_bytes()).hexdigest(),"test":hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}}
    (HERE/"plan/local-two-job-smoke.json").write_text(json.dumps(result,indent=2,sort_keys=True)+"\n");print(json.dumps(result,sort_keys=True))
if __name__=="__main__":main()
