#!/usr/bin/env python3
"""Path-parameterized exact audit of the archived n=6 h1..h5 residues."""
from __future__ import annotations
import argparse, hashlib, json, math, re
from fractions import Fraction
from pathlib import Path

REPS = {0: 4, 1: 8, 2: 8, 7: 4, 8: 8, 14: 4}
FIELDS = [f"h{i}_residue" for i in range(1, 6)]

def sha(p: Path) -> str: return hashlib.sha256(p.read_bytes()).hexdigest()
def crt(a, ps):
    q=math.prod(ps); return sum(x*(q//p)*pow(q//p,-1,p) for x,p in zip(a,ps))%q
def garner(a, ps):
    x, b=0,1
    for r,p in zip(a,ps): x += b*(((r-x)%p)*pow(b%p,-1,p)%p); b *= p
    return x
def pratt(nodes, k, active):
    p=int(k)
    if p==2:return True
    if p<3 or p%2==0 or k in active or k not in nodes:return False
    r=nodes[k]; fac={int(q):int(e) for q,e in r['p_minus_1_factorization'].items()}
    if math.prod(q**e for q,e in fac.items())!=p-1:return False
    a=int(r['witness'])
    if not 1<a<p or pow(a,p-1,p)!=1:return False
    active.add(k); ok=all(pratt(nodes,str(q),active) and math.gcd(pow(a,(p-1)//q,p)-1,p)==1 for q in fac); active.remove(k)
    return ok
def resource(p:Path):
    t=p.read_text(); m=re.search(r"^\s*([0-9.]+) real\s+([0-9.]+) user\s+([0-9.]+) sys$",t,re.M)
    rss=re.search(r"^\s*(\d+)  maximum resident set size$",t,re.M); swaps=re.search(r"^\s*(\d+)  swaps$",t,re.M)
    assert m and rss and swaps
    return {'real_seconds':m.group(1),'user_seconds':m.group(2),'sys_seconds':m.group(3),'rss_bytes':int(rss.group(1)),'swaps':int(swaps.group(1))}
def main():
    p=argparse.ArgumentParser(); p.add_argument('--runs',type=Path,required=True); p.add_argument('--source',type=Path,required=True); p.add_argument('--pratt',type=Path,required=True); p.add_argument('--punctured',type=Path,required=True); p.add_argument('--top2',type=Path,required=True); p.add_argument('--output',type=Path,required=True); a=p.parse_args()
    cert=json.loads(a.pratt.read_text()); ps=list(map(int,cert['roots'])); assert len(ps)==4 and all(pratt(cert['nodes'],str(x),set()) for x in ps)
    bound=math.factorial(35); assert math.prod(ps)>bound
    reconstructed={}; records={}; peak=0
    for v,mult in REPS.items():
        rows=[]; rec=[]
        for i,prime in enumerate(ps):
            stem=a.runs/f'n6-v{v}-p{i}'; x=json.loads(stem.with_suffix('.json').read_text()); rr=resource(stem.with_suffix('.resource.txt'))
            assert stem.with_suffix('.exit-status.txt').read_text()=='exit_status=0\n'
            assert x['n']==6 and x['v']==v and x['prime_index']==i and int(x['prime'])==prime and x['threads']==12
            assert x['algorithm']=='short-gap-h1-to-h5-modular-layer-rolling' and rr['swaps']==0
            rows.append(x); peak=max(peak,rr['rss_bytes']); rec.append({'json':stem.with_suffix('.json').name,'json_sha256':sha(stem.with_suffix('.json')),'exit':stem.with_suffix('.exit-status.txt').name,'exit_sha256':sha(stem.with_suffix('.exit-status.txt')),'resource':stem.with_suffix('.resource.txt').name,'resource_sha256':sha(stem.with_suffix('.resource.txt')),'summary':rr})
        h={}
        for f in FIELDS:
            rs=[int(x[f]) for x in rows]; x,y=crt(rs,ps),garner(rs,ps); assert x==y and 0<=x<=bound and all(x%q==z for q,z in zip(ps,rs)); h[f.removesuffix('_residue')]=x
        reconstructed[v]=h; records[str(v)]={'multiplicity':mult,'h':h,'files':rec}
    h={f'h{i}':sum(REPS[v]*reconstructed[v][f'h{i}'] for v in REPS) for i in range(1,6)}
    punct=json.loads(a.punctured.read_text()); top=json.loads(a.top2.read_text()); pv={int(x['v']):x for x in punct['per_vertex']}
    B=sum(REPS[v]*int(pv[v]['B_v']) for v in REPS); A2=sum(REPS[v]*int(pv[v]['A2_v']) for v in REPS); A1=int(top['A1'])
    assert B==int(punct['B'])==int(top['Btop2']) and A2==int(punct['A2'])==int(top['A2'])
    assert 4*h['h1']-2*B==A1 and h['h2']==2*h['h3']
    tails={}; subtotal=0
    for t in range(2,7):
        subtotal+=h[f'h{t-1}']; z=B-subtotal; frac=Fraction(z,B); assert z>=0
        tails[f'T{t}']={'count':z,'fraction_of_B':{'numerator':frac.numerator,'denominator':frac.denominator},'QE8_holds':t*(t+1)*z<=8*B}
    out={'status':'passed','scope':'strict finite n=6 CRT reconstruction; no asymptotic inference','source_sha256':sha(a.source),'pratt_sha256':sha(a.pratt),'references':{'punctured_sha256':sha(a.punctured),'top2_sha256':sha(a.top2)},'parameters':{'n':6,'representatives':REPS,'threads':12,'primes':ps,'modulus':math.prod(ps),'bound':bound},'runs':records,'aggregate_h':h,'B':B,'A1':A1,'A2':A2,'reference_B':B,'reference_A1':A1,'reference_A2':A2,'checks':{'A1_equals_4h1_minus_2B':True,'h2_equals_2h3':True,'all_swaps_zero':True,'peak_rss_bytes':peak},'tails':tails}
    a.output.parent.mkdir(parents=True,exist_ok=True); a.output.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n'); print(json.dumps({'status':'passed','aggregate_h':h,'peak_rss_bytes':peak},indent=2))
if __name__=='__main__':main()
