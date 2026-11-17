#!/usr/bin/env python3
"""Verify 40 residues, D4 orbits, and independent global integer totals."""
import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import subprocess
import sys

ap=argparse.ArgumentParser()
ap.add_argument('--remote-root',type=Path,required=True)
ap.add_argument('--r028-layer',type=Path,required=True)
ap.add_argument('--global-reference',type=Path,required=True)
ap.add_argument('--outdir',type=Path,required=True)
args=ap.parse_args()
args.outdir.mkdir(exist_ok=True)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
source_hash='5f93fa39695e0f4badabcbabaad2befab920fa1590ea2337fe47a09201587db9'
binary_hash='e57308e1d422b949ed43c39f8939f3fb046d7041e6f1adf61a93d48b7994c37c'
for part in ('v0','followon'):
    assert sha(args.remote_root/part/'punctured_layer_mod.cpp')==source_hash
    assert sha(args.remote_root/part/'punctured_layer_mod')==binary_hash
assert sha(args.global_reference)=='0d3ac2e832e4fb37b14d5ed157ef8036811b88be7123410fa465781c9bbfdbf6'
ref=json.loads(args.global_reference.read_text())
cert=args.r028_layer/'n7-center-crt/pratt-certificates.json'
verifier=Path(__file__).with_name('verify_selected_hole.py')
seen=set(); table=[]; exit_checks=[]
for v in (0,7,8,14,15,16,21,22,23,24):
    x,y=v%7,v//7
    orbit={b*7+a for a,b in ((x,y),(y,x),(6-x,y),(x,6-y),(6-x,6-y),(6-y,x),(y,6-x),(6-y,6-x))}
    assert not seen & orbit
    seen |= orbit
    if v==24:
        paths=[args.r028_layer/f'n7-center-crt/residues/n7-v24-p{i}.json' for i in range(4)]
    elif v==0:
        paths=[args.remote_root/f'v0/n7-v0-p{i}.json' for i in range(4)]
    else:
        paths=[args.remote_root/f'followon/results/v{v}/p{i}.json' for i in range(4)]
    if v!=24:
        for i,p in enumerate(paths):
            receipt=p.with_suffix('.exit.json')
            resource=p.with_suffix('.resource.txt')
            r=json.loads(receipt.read_text())
            assert r['exit_code']==0 and r['output_sha256']==sha(p)
            assert r['resource_sha256']==sha(resource)
            exit_checks.append({'v':v,'prime_index':i,'output_sha256':sha(p),'receipt_sha256':sha(receipt),'resource_sha256':sha(resource),'exit_code':0})
    result=subprocess.run([sys.executable,'-B',str(verifier),'--vertex',str(v),str(cert),*map(str,paths)],capture_output=True,text=True,check=True)
    d=json.loads(result.stdout)
    assert d['verified'] and d['v']==v
    (args.outdir/f'v{v}-verified.json').write_text(result.stdout)
    A,B=int(d['A2v']),int(d['Bv'])
    ratio=Fraction(B,A)
    table.append({'v':v,'xy':[x,y],'orbit':sorted(orbit),'multiplicity':len(orbit),'A2v':str(A),'Bv':str(B),'Bv_over_A2v_exact':str(ratio),'Bv_over_A2v_display':float(ratio),'above_one_fifth':5*B>=A})
assert seen==set(range(49))
A=sum(r['multiplicity']*int(r['A2v']) for r in table)
B=sum(r['multiplicity']*int(r['Bv']) for r in table)
assert A==int(ref['A2']) and B==int(ref['Btop2'])
lowest=min(table,key=lambda r:Fraction(int(r['Bv']),int(r['A2v'])))
out={'verified':True,'n':7,'scope':'exact finite full D4 punctured table; uniform-in-n bounds remain unproved','source_sha256':source_hash,'binary_sha256':binary_hash,'global_reference_sha256':sha(args.global_reference),'pratt_certificate_sha256':sha(cert),'vertex_verifier_sha256':sha(verifier),'orchestrator_sha256':sha(Path(__file__)),'table':table,'sum_A2v':str(A),'sum_Bv':str(B),'global_B_over_A2_display':B/A,'global_totals_match_independent_reference':True,'minimum_ratio_vertex':lowest['v'],'minimum_ratio_display':lowest['Bv_over_A2v_display'],'all_vertices_at_least_one_fifth':all(r['above_one_fifth'] for r in table),'new_pass_exit_checks':exit_checks,'prior_center':'r028 exact residues reused, not rerun'}
(args.outdir/'n7-orbits-verified.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:v for k,v in out.items() if k not in ('table','new_pass_exit_checks')},indent=2))
