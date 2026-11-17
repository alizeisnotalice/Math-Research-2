#!/usr/bin/env python3
"""Re-run the archived exact identity checks and both squarefree certificates."""
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import hashlib, json, platform, shutil, subprocess, sys, time
BASE=Path(__file__).resolve().parent
jobs=[(BASE/"results/041-slice"/n,[sys.executable]) for n in ["audit_quartic_factor_types.py","check_analytic_two_by_two_closure.py","check_reducible_quartic_terminal.py"]]
jobs.append((BASE/"results/041-squarefree/verify_q2_fraction_interval.py",[sys.executable]))
wl=shutil.which("wolframscript")
if wl is None: raise SystemExit("wolframscript unavailable: no Wolfram verification claimed")
jobs.append((BASE/"results/041-squarefree/verify_q2_newman_squarefree.wl",[wl,"-file"]))
def run(job):
    path,cmd=job;start=time.monotonic();r=subprocess.run([*cmd,str(path)],capture_output=True,text=True)
    return dict(file=str(path.relative_to(BASE)),source_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),returncode=r.returncode,seconds=time.monotonic()-start,stdout=r.stdout,stderr=r.stderr)
with ThreadPoolExecutor(max_workers=5) as pool: rows=list(pool.map(run,jobs))
version=subprocess.run([wl,"-code","$Version"],capture_output=True,text=True)
report=dict(all_passed=all(r["returncode"]==0 for r in rows) and version.returncode==0,python=sys.version,platform=platform.platform(),wolfram_version=version.stdout.strip(),checks=rows,scope="Three exact symbolic regressions and two independently implemented finite-state certificates; original conjecture remains open.")
(BASE/"results/041-proof-parent-checks.json").write_text(json.dumps(report,indent=2,ensure_ascii=False)+"\n")
print(json.dumps({"all_passed":report["all_passed"],"wolfram_version":report["wolfram_version"],"checks":[{k:r[k] for k in ("file","returncode","seconds")} for r in rows]},ensure_ascii=False))
raise SystemExit(0 if report["all_passed"] else 1)
