#!/usr/bin/env python3
"""Path-parameterized finite n=7 centre CRT audit; no production rerun."""
from __future__ import annotations
import argparse, hashlib, json, math, re
from pathlib import Path
FIELDS=[f'h{i}_residue' for i in range(1,6)]
STEMS=['retry1-n7-v24-p0','n7-v24-p1','n7-v24-p2','n7-v24-p3']
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def pratt(nodes,k,active):
 p=int(k)
 if p==2:return True
 if p<3 or p%2==0 or k in active or k not in nodes:return False
 r=nodes[k]; f={int(q):int(e) for q,e in r['p_minus_1_factorization'].items()}
 if math.prod(q**e for q,e in f.items())!=p-1:return False
 a=int(r['witness'])
 if not 1<a<p or pow(a,p-1,p)!=1:return False
 active.add(k); ok=all(pratt(nodes,str(q),active) and math.gcd(pow(a,(p-1)//q,p)-1,p)==1 for q in f); active.remove(k); return ok
def crt(a,ps):
 q=math.prod(ps); return sum(x*(q//p)*pow(q//p,-1,p) for x,p in zip(a,ps))%q
def garner(a,ps):
 x,b=0,1
 for r,p in zip(a,ps):x+=b*(((r-x)%p)*pow(b%p,-1,p)%p); b*=p
 return x
def resource(p):
 t=Path(p).read_text(); m=re.search(r"^\s*([0-9.]+) real\s+([0-9.]+) user\s+([0-9.]+) sys$",t,re.M); rss=re.search(r"^\s*(\d+)  maximum resident set size$",t,re.M); swaps=re.search(r"^\s*(\d+)  swaps$",t,re.M); assert m and rss and swaps
 return {'real_seconds':m.group(1),'user_seconds':m.group(2),'sys_seconds':m.group(3),'rss_bytes':int(rss.group(1)),'swaps':int(swaps.group(1))}
def main():
 p=argparse.ArgumentParser(); p.add_argument('--runs',type=Path,required=True); p.add_argument('--source',type=Path,required=True); p.add_argument('--pratt',type=Path,required=True); p.add_argument('--old-center',type=Path,required=True); p.add_argument('--normalizer',type=Path,required=True); p.add_argument('--output',type=Path,required=True); a=p.parse_args()
 cert=json.loads(a.pratt.read_text()); ps=list(map(int,cert['roots'])); assert len(ps)==4 and all(pratt(cert['nodes'],str(x),set()) for x in ps); bound=math.factorial(48); assert math.prod(ps)>bound
 rows=[]; files=[]; peak=0
 for i,stem0 in enumerate(STEMS):
  stem=a.runs/stem0; x=json.loads(stem.with_suffix('.json').read_text()); rr=resource(stem.with_suffix('.resource.txt')); assert stem.with_suffix('.exit-status.txt').read_text()=='exit_status=0\n'; assert x['n']==7 and x['v']==24 and x['prime_index']==i and int(x['prime'])==ps[i] and x['threads']==12 and x['algorithm']=='short-gap-h1-to-h5-modular-layer-rolling' and rr['swaps']==0
  rows.append(x); peak=max(peak,rr['rss_bytes']); files.append({'prime_index':i,'stem':stem0,'json_sha256':sha(stem.with_suffix('.json')),'exit_sha256':sha(stem.with_suffix('.exit-status.txt')),'resource_sha256':sha(stem.with_suffix('.resource.txt')),'summary':rr})
 h={}
 for f in FIELDS:
  rs=[int(x[f]) for x in rows]; x,y=crt(rs,ps),garner(rs,ps); assert x==y and 0<=x<=bound and all(x%q==z for q,z in zip(ps,rs)); h[f.removesuffix('_residue')]=x
 old=json.loads(a.old_center.read_text()); norm=json.loads(a.normalizer.read_text()); assert all(h[f'h{i}']==int(old['h'][f'h{i}']) for i in range(1,4)); assert int(old['reference_selected_normalizer']['Bv'])==int(norm['Bv']) and int(old['reference_selected_normalizer']['A2v'])==int(norm['A2v'])
 relation=h['h3']-3*h['h4']+2*h['h5']
 out={'status':'passed','scope':'strict finite n=7 centre CRT reconstruction only; no asymptotic inference','source_sha256':sha(a.source),'pratt_sha256':sha(a.pratt),'references':{'old_center_sha256':sha(a.old_center),'normalizer_sha256':sha(a.normalizer)},'parameters':{'n':7,'v':24,'threads':12,'parallel_processes':1,'randomness':'none','primes':ps,'modulus':math.prod(ps),'bound':bound},'h':h,'raw_files':files,'reference':{'Bv':int(norm['Bv']),'A2v':int(norm['A2v']),'old_h1_h3':{f'h{i}':int(old['h'][f'h{i}']) for i in range(1,4)}},'checks':{'h1_h2_h3_equal_old_center':True,'selected_normalizer_Bv_A2v_checked':True,'all_swaps_zero':True,'peak_rss_bytes':peak},'record_only_relation':{'expression':'h3 - 3*h4 + 2*h5','value':relation,'note':'recorded exactly; no per-hole identity is asserted'}}
 a.output.parent.mkdir(parents=True,exist_ok=True); a.output.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n'); print(json.dumps({'status':'passed','h':h,'relation':relation,'peak_rss_bytes':peak},indent=2))
if __name__=='__main__':main()
