#!/usr/bin/env python3
"""Portable guard for the current calibration-derived source smoke trace."""
import hashlib,json,sys
from pathlib import Path
SOURCE='3a9f3c0a8c10409346c25f16f2d82463f928bfaacf7307b4a88a770ec8f69ec6'
def main():
 raw,trace,src=map(Path,sys.argv[1:4]); assert hashlib.sha256(src.read_bytes()).hexdigest()==SOURCE
 lines=[json.loads(x) for x in raw.read_text().splitlines()]; meta=lines[0]; assert meta['implementation']=='r026-tau2-firstmerge-rb-terminal'
 rb=next(x for x in lines if x['kind']=='firstmerge_rb_calibration'); t=json.loads(trace.read_text()); assert t['n']==3 and t['P']==17
 p=next(x for x in t['rows'] if x.get('terminal'))['first_merge_rb_payload_ab_r1_r2']; n=t['n']; P=t['P']; out=[0.0]*4
 for i in range(P):
  a,b,x,y=p[4*i:4*i+4];q=a+b;m=q*(q-1)/(2*a*b);dx=x%n-y%n;dy=x//n-y//n;d2=dx*dx+dy*dy;ax=2*(x%n)-(n-1);ay=2*(x//n)-(n-1);bx=2*(y%n)-(n-1);by=2*(y//n)-(n-1);f=ax*ax+ay*ay+bx*bx+by*by;g=ax*bx+ay*by
  assert q*(q-1)*(abs(dx)+abs(dy))>=4*a*b*(abs(dx)+abs(dy)-1) and f-2*g==4*d2
  for j,v in enumerate((m,m*d2,m*f,m*g)):out[j]+=v/P
 for k,v in zip(('mean_m','mean_m_d2','mean_m_F0','mean_m_G0'),out):assert abs(rb[k]-v)<1e-12
 print(json.dumps({'audit_passed':True,'scope':'derived source hash and trace GPU payload reproduce new m fields; not an accuracy claim'}))
if __name__=='__main__':main()
