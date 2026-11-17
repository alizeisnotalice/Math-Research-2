#!/usr/bin/env python3
"""Exact small audits; does not launch production search or graph scans."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import time

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--out',type=Path,required=True)
    a=p.parse_args()
    root=Path(__file__).resolve().parent
    out=a.out.resolve()
    out.mkdir(parents=True,exist_ok=False)
    jobs=[('ordered-branch','audit_n3_ordered_kernel.py',[]),
          ('singleton-geometry','audit_geometric_psi_n3.py',[]),
          ('singleton-geometry','audit_far_singleton_family.py',[]),
          ('singleton-geometry','audit_high_singleton_identity_n3.py',[]),
          ('early-low-branch','audit_corner_swap_n3.py',[]),
          ('binomial-crosscheck','audit_fibres.py',[]),
          ('finite-candidates','verify_psi.py',['n8-bfs-center-s2.txt']),
          ('finite-candidates','verify_psi.py',['n12-bfs-center-s1.txt'])]
    for group in sorted({j[0] for j in jobs}):
        shutil.copytree(root/group,out/group)
    records=[]
    for i,(group,name,args) in enumerate(jobs):
        source=out/group/name
        began=time.time()
        with (out/f'check-{i}.stdout').open('w') as f, (out/f'check-{i}.stderr').open('w') as e:
            child=subprocess.Popen([sys.executable,'-B',str(source),*args],
                                   cwd=out/group,stdout=f,stderr=e)
            rc=child.wait()
        absent=subprocess.run(['ps','-p',str(child.pid),'-o','pid='],capture_output=True).returncode!=0
        records.append(dict(group=group,script=name,args=args,
                            source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
                            pid=child.pid,returncode=rc,waited_reaped=True,
                            absent_after=absent,seconds=time.time()-began))
        (out/'receipt.json').write_text(json.dumps(dict(python=sys.version,checks=records),indent=2)+'\n')
        assert rc==0 and absent,(group,name,rc)
    print(json.dumps(dict(status='passed',checks=len(records),production_rerun=False)))

if __name__=='__main__':
    main()
