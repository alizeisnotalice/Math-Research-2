#!/usr/bin/env python3
import copy
import hashlib
import importlib.util
import json
import math
import tempfile
import os
from pathlib import Path

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location("audit",HERE/"aggregate_production_r026.py")
audit=importlib.util.module_from_spec(spec);spec.loader.exec_module(audit)
CERT_ROOT=Path(os.environ.get('A022_CERTIFICATES_ROOT',HERE.parent.parent))
CAL=CERT_ROOT/'r026'/'r026-merge-rb-calibration'
source=CAL/"source/tau2_merge_rb_calibration_r026.swift"
results={}
for n,filename in [(3,"n3-p4096-r64.jsonl"),(7,"n7-p16384-r64.jsonl")]:
    result=audit.audit(CAL/"raw"/filename,source)
    old=json.loads((CAL/f"results/n{n}-summary.json").read_text())
    assert math.isclose(result["Q_mRB"],old["Q_mRB"],rel_tol=2e-14)
    assert math.isclose(result["Q_L_legacy"],old["Q_L_legacy"],rel_tol=2e-14)
    (HERE/f"n{n}-production-audit-final.json").write_text(json.dumps(result,indent=2)+"\n")
    results[str(n)]={"Q_mRB":result["Q_mRB"],"Q_L_legacy":result["Q_L_legacy"],"raw_sha256":result["raw_sha256"]}
original=[json.loads(s) for s in (CAL/"raw/n3-p4096-r64.jsonl").read_text().splitlines()]
with tempfile.TemporaryDirectory(prefix="audit-cases-",dir=HERE) as td:
    path=Path(td)/"case.jsonl"
    def evaluate(rows):
        path.write_text("\n".join(json.dumps(x) for x in rows)+"\n")
        return audit.audit(path,source)
    rejected={}
    def reject(name,rows):
        try:evaluate(rows)
        except (ValueError,KeyError,TypeError) as exc:rejected[name]=str(exc)
        else:raise AssertionError("accepted "+name)
    a=copy.deepcopy(original);a[4]=copy.deepcopy(a[2]);reject("duplicate RB id",a)
    a=copy.deepcopy(original);a.pop(2);reject("missing RB",a)
    a=copy.deepcopy(original);a[0]["implementation"]="wrong";reject("wrong implementation",a)
    a=copy.deepcopy(original);a[1]["seed_hex"]="0";reject("wrong seed",a)
    a=copy.deepcopy(original);a[1]["logZ"]=float("nan");reject("nonfinite logZ",a)
    a=copy.deepcopy(original);a[1]["Z"]=float("inf");reject("infinite Z",a)
    a=copy.deepcopy(original);a[1]["tau2_lift_L_G0_groups"][0]+=1;reject("integer identity broken",a)
    a=copy.deepcopy(original);a[2]["mean_m_hist"][2]+=.1;reject("RB histogram broken",a)
    reject("missing run_end",copy.deepcopy(original[:-1]))
    reordered=evaluate([original[0],*reversed(original[1:-1]),original[-1]])
    assert reordered["Q_mRB"]==results["3"]["Q_mRB"]
    a=copy.deepcopy(original)
    for row in a:
        if row.get("kind")=="replicate":row.update(logZ=-1000.0-row["replicate"],Z=0.0)
    underflow=evaluate(a)
    assert underflow["B_over_Nfactorial_display"]==0 and math.isfinite(underflow["log_B_over_Nfactorial"])
    results.update(rejected=rejected,reordered_records_accepted=True,synthetic_underflow={"scope":"aggregator fixture only, not an n32 experiment","log_B_over_Nfactorial":underflow["log_B_over_Nfactorial"],"display":underflow["B_over_Nfactorial_display"]})
results["auditor_sha256"]=hashlib.sha256((HERE/"aggregate_production_r026.py").read_bytes()).hexdigest()
(HERE/"production-audit-tests.json").write_text(json.dumps(results,indent=2)+"\n")
print(json.dumps({"passed":True,"negative_cases":len(rejected),"auditor_sha256":results["auditor_sha256"]}))
