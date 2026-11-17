#!/usr/bin/env python3
"""Re-run the three archived exact checks and record exits and source hashes."""
import hashlib, json, platform, subprocess, sys, time
from pathlib import Path
HERE=Path(__file__).resolve().parent
rows=[]
for name in ["check_final_low_degree_slice_closure.py", "counterexample_exact.py", "offp_overlap_verify_phi6.py"]:
    start=time.monotonic()
    run=subprocess.run([sys.executable,str(HERE/name)],capture_output=True,text=True)
    rows.append({"script":name,"source_sha256":hashlib.sha256((HERE/name).read_bytes()).hexdigest(),"returncode":run.returncode,"seconds":time.monotonic()-start,"stdout":run.stdout,"stderr":run.stderr})
report={"python":sys.version,"platform":platform.platform(),"checks":rows,"all_passed":all(r["returncode"]==0 for r in rows),"scope":"Exact regression identities and two explicit local-lemma counterexamples; not a proof of the original conjecture."}
(HERE/"parent-recheck.json").write_text(json.dumps(report,indent=2,ensure_ascii=False)+"\n")
print(json.dumps({"all_passed":report["all_passed"],"checks":[{k:r[k] for k in ("script","returncode","seconds")} for r in rows]}))
raise SystemExit(0 if report["all_passed"] else 1)
