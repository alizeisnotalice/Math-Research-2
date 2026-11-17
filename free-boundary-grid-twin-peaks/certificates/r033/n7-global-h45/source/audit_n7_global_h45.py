#!/usr/bin/env python3
"""Independent exact reconstruction of all ten n=7 D4 h1..h5 rows."""
from __future__ import annotations
import argparse, hashlib, json, math, re
from fractions import Fraction
from pathlib import Path

MULT={0:4,7:8,8:4,14:8,15:8,16:4,21:4,22:4,23:4,24:1}
FIELDS=[f'h{i}_residue' for i in range(1,6)]
def sha(p:Path)->str:return hashlib.sha256(p.read_bytes()).hexdigest()
def pratt(nodes,k,active):
 p=int(k)
 if p==2:return True
 if p<3 or p%2==0 or k in active or k not in nodes:return False
 r=nodes[k]; fac={int(q):int(e) for q,e in r['p_minus_1_factorization'].items()}
 if math.prod(q**e for q,e in fac.items())!=p-1:return False
 a=int(r['witness'])
 if not 1<a<p or pow(a,p-1,p)!=1:return False
 active.add(k); ok=all(pratt(nodes,str(q),active) and math.gcd(pow(a,(p-1)//q,p)-1,p)==1 for q in fac); active.remove(k); return ok
def crt(a,ps):
 q=math.prod(ps); return sum(x*(q//p)*pow(q//p,-1,p) for x,p in zip(a,ps))%q
def garner(a,ps):
 x,b=0,1
 for r,p in zip(a,ps):x+=b*(((r-x)%p)*pow(b%p,-1,p)%p);b*=p
 return x
def resource(p:Path):
 t=p.read_text(); rss=re.search(r'^\s*(\d+)  maximum resident set size$',t,re.M); swaps=re.search(r'^\s*(\d+)  swaps$',t,re.M); tm=re.search(r'^\s*([0-9.]+) real\s+([0-9.]+) user\s+([0-9.]+) sys$',t,re.M); assert rss and swaps and tm
 return {'rss_bytes':int(rss.group(1)),'swaps':int(swaps.group(1)),'real_seconds':tm.group(1),'user_seconds':tm.group(2),'sys_seconds':tm.group(3)}
def fraction(x,y):
 q=Fraction(x,y); return {'numerator':q.numerator,'denominator':q.denominator,'decimal_display_only':float(q)}
def main():
 p=argparse.ArgumentParser();p.add_argument('--raw',type=Path,required=True);p.add_argument('--source',type=Path,required=True);p.add_argument('--pratt',type=Path,required=True);p.add_argument('--r31',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
 cert=json.loads(a.pratt.read_text()); ps=list(map(int,cert['roots'])); assert len(ps)==4 and all(pratt(cert['nodes'],str(x),set()) for x in ps); bound=math.factorial(48); assert math.prod(ps)>bound
 r31=json.loads(a.r31.read_text()); assert r31['status']=='passed'; ref_orbits=r31['orbits']; ref_global=r31['global']; assert sum(MULT.values())==49
 out_orbits={}; peak=0
 for v,m in MULT.items():
  rows=[]; files=[]
  for i,prime in enumerate(ps):
   stem=a.raw/f'v{v}'/f'n7-v{v}-p{i}'; j=stem.with_suffix('.json'); x=json.loads(j.read_text()); rp=stem.with_suffix('.resource.txt'); rr=resource(rp); assert x['n']==7 and x['v']==v and x['prime_index']==i and int(x['prime'])==prime and x['algorithm']=='short-gap-h1-to-h5-modular-layer-rolling' and rr['swaps']==0
   if v in (22,23):
    ep=stem.with_suffix('.exit.json'); receipt=json.loads(ep.read_text()); assert receipt['returncode']==0 and receipt['output_sha256']==sha(j) and receipt['resource_sha256']==sha(rp)
   else:
    ep=stem.with_suffix('.exit-status.txt'); assert ep.read_text()=='exit_status=0\n'
   peak=max(peak,rr['rss_bytes']); rows.append(x); files.append({'prime_index':i,'json_sha256':sha(j),'exit_sha256':sha(ep),'resource_sha256':sha(rp),'resource':rr})
  h={}
  for f in FIELDS:
   residues=[int(x[f]) for x in rows]; x,y=crt(residues,ps),garner(residues,ps); assert x==y and 0<=x<=bound and all(x%q==z for q,z in zip(ps,residues));h[f.removesuffix('_residue')]=x
  old=ref_orbits[str(v)]; assert int(old['multiplicity'])==m and all(h[f'h{i}']==int(old['h'][f'h{i}']) for i in range(1,4))
  out_orbits[str(v)]={'multiplicity':m,'h':h,'files':files,'r31_h1_h3_match':True}
 h={f'h{i}':sum(MULT[v]*out_orbits[str(v)]['h'][f'h{i}'] for v in MULT) for i in range(1,6)}
 B,A2=int(ref_global['B']),int(ref_global['A2']); assert all(h[f'h{i}']==int(ref_global['h'][f'h{i}']) for i in range(1,4))
 tails={}; partial=0
 for t in range(2,7):
  partial+=h[f'h{t-1}']; z=B-partial; assert z>=0; tails[f'T{t}']={'count':z,'B_fraction':fraction(z,B),'QE8_lhs':t*(t+1)*z,'QE8_rhs':8*B,'QE8_slack':8*B-t*(t+1)*z,'QE8_holds':t*(t+1)*z<=8*B}
 unseenB=B-sum(h.values()); unseenA=A2-sum(i*h[f'h{i}'] for i in range(1,6)); assert unseenB>0 and unseenA>=6*unseenB
 result={'status':'passed','scope':'exact finite n=7 all-D4-orbit h1..h5 reconstruction; no asymptotic conclusion','source_sha256':sha(a.source),'pratt_sha256':sha(a.pratt),'r31_reference_sha256':sha(a.r31),'parameters':{'n':7,'orbits':MULT,'multiplicity_sum':sum(MULT.values()),'primes':ps,'crt_modulus':math.prod(ps),'count_bound':bound},'orbits':out_orbits,'global':{'B':B,'A2':A2,'h':h,'r31_h1_h3_match':True,'h3_minus_3h4_plus_2h5':h['h3']-3*h['h4']+2*h['h5'],'tails_T2_to_T6':tails,'L_at_least_6':{'B_fraction':fraction(unseenB,B),'lifted_fraction':fraction(unseenA,A2),'conditional_mean':fraction(unseenA,unseenB)}},'checks':{'all_exit_zero':True,'all_swaps_zero':True,'peak_rss_bytes':peak,'D4_multiplicity_sum_49':True}}
 a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n');print(json.dumps({'status':'passed','h':h,'R45':result['global']['h3_minus_3h4_plus_2h5'],'T6':tails['T6'],'L_ge_6_lifted':result['global']['L_at_least_6']['lifted_fraction']},indent=2))
if __name__=='__main__':main()
