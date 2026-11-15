#!/usr/bin/env python3
"""Exact defining equation/positive branch check for the actual R52 raw wall."""
from fractions import Fraction as F
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import platform
import sys
import time

sys.dont_write_bytecode=True
begin=time.process_time();inputpath=Path(sys.argv[1]);witnesspath=Path(sys.argv[2]);outpath=Path(sys.argv[3]);raw=inputpath.read_bytes();d=json.loads(raw);w=json.loads(witnesspath.read_bytes())
module=Path(__file__).with_name('r057_quadratic_endpoint_kernel.py');spec=importlib.util.spec_from_file_location('r057kernel',module);k=importlib.util.module_from_spec(spec);spec.loader.exec_module(k)
field=k.Field(d['Q']);event=field.element(*d['event']);t0=F(d['t0']);A,B=[r['expr'] for r in d['endpoints']]
assert int(A[2])==-1 and int(B[2])==1 and F(A[4])==0
Q=F(A[3]);assert Q==field.Q
D=F(A[0])-F(B[0]);q0=F(B[3])-F(B[4])*t0;q1=F(B[4]);assert q1<0
alpha,beta=(D*D+Q-q0)/q1,-2*D/q1
assert event==field.element(alpha,beta)
branch=field.element(D,-1);assert branch.sign()>0
UB=field.element(q0)+q1*event;assert UB==branch*branch and UB.sign()>0
assert event.sign()>0 and (F(1,2)-event).sign()>0
assert (event-F(w['t_lower'])).sign()>0 and (F(w['t_upper'])-event).sign()>0
left,right=k.prepare_recipe(A,t0,event),k.prepare_recipe(B,t0,event);relation=k.contact(left,right)
assert relation['center_sign']==0 and relation['contact_order']==1 and relation['left_sign']==-1 and relation['right_sign']==1
out={'status':'passed','scope':'actual two-endpoint event equation, positive root branch, exact event enclosure and first jet; no exposed-union-wall assertion',
 'event':event.pair(),'Q':str(Q),'positive_branch_D_minus_sqrtQ':branch.pair(),'radicand_B_at_event':UB.pair(),'relation_A_minus_B':relation,
 'input_sha256':hashlib.sha256(raw).hexdigest(),'witness_sha256':hashlib.sha256(witnesspath.read_bytes()).hexdigest(),
 'kernel_sha256':hashlib.sha256(module.read_bytes()).hexdigest(),'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
 'python':sys.version,'platform':platform.platform(),'cpu_seconds':time.process_time()-begin,'actual_command':[sys.executable,*sys.argv]}
if not os.path.ismount('/Volumes/data'):raise RuntimeError('/Volumes/data is not mounted')
outpath.parent.mkdir(parents=True,exist_ok=True);outpath.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
