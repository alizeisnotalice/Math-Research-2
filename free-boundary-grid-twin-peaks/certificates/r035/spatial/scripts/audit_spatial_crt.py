#!/usr/bin/env python3
"""Exact finite CRT audit of the 40 n=7 spatial-kernel residues."""
from __future__ import annotations
import argparse, hashlib, json, math
from fractions import Fraction
from pathlib import Path

MULT={0:4,7:8,8:4,14:8,15:8,16:4,21:4,22:4,23:4,24:1}
PRIMES=[9223372036854775783,9223372036854775643,9223372036854775549,9223372036854775507]
FIELDS=[f'h{i}_residue' for i in range(1,6)]+[f'E{i}_residue' for i in range(1,6)]+[f'far{k}_h{i}_residue' for k in range(1,4) for i in range(1,6)]
SOURCE='e1bc6ec7753f60f65f61da7b332e481c8cbe4613b1a31f6147935da081a52595'
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def pratt(nodes,k,active):
 p=int(k)
 if p==2:return True
 if p<3 or p%2==0 or k in active or k not in nodes:return False
 r=nodes[k]; fac={int(q):int(e) for q,e in r['p_minus_1_factorization'].items()}
 if math.prod(q**e for q,e in fac.items())!=p-1:return False
 a=int(r['witness'])
 if not 1<a<p or pow(a,p-1,p)!=1:return False
 active.add(k); ok=all(pratt(nodes,str(q),active) and math.gcd(pow(a,(p-1)//q,p)-1,p)==1 for q in fac);active.remove(k);return ok
def crt(a,ps):
 q=math.prod(ps);return sum(x*(q//p)*pow(q//p,-1,p) for x,p in zip(a,ps))%q
def garner(a,ps):
 x,b=0,1
 for r,p in zip(a,ps): x+=b*(((r-x)%p)*pow(b%p,-1,p)%p);b*=p
 return x
def frac(x,y):
 z=Fraction(x,y);return {'numerator':z.numerator,'denominator':z.denominator,'decimal_display_only':float(z)}
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[1]);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();root=a.root.resolve()
 assert sha(root/'source/gap_spatial_kernel_mod.cpp')==SOURCE
 cert=json.loads((root/'references/pratt-certificates.json').read_text());assert list(map(int,cert['roots']))==PRIMES and all(pratt(cert['nodes'],str(p),set()) for p in PRIMES)
 fact=math.factorial(48);mod=math.prod(PRIMES);assert mod>72*fact and sum(MULT.values())==49
 source_audit=json.loads((root/'evidence/source-complete-audit.json').read_text());assert source_audit['source_sha256']==SOURCE and len(source_audit['runs'])==40 and source_audit['checks']['all_rc0'] and source_audit['checks']['all_scalar_h_reference_checked'] and source_audit['checks']['all_swaps_zero']
 old=json.loads((root/'references/r33-n7-global-h45-crt-audit.json').read_text());assert old['status']=='passed'
 out={};glob={f[:-8]:0 for f in FIELDS}; nonneg={};support={};r45={}
 for v,m in MULT.items():
  rows=[]
  for i,p in enumerate(PRIMES):
   x=json.loads((root/f'raw/n7-v{v}-p{i}.json').read_text());assert (x['n'],x['v'],x['prime_index'],int(x['prime']))==(7,v,i,p) and x['algorithm']=='gap-spatial-h1-to-h5-modular-layer-rolling' and x['n7_required_explicit_optin'] and x['max_modeled_payload_bytes']<=x['payload_limit_bytes'];rows.append(x)
  vals={}
  for f in FIELDS:
   rr=[int(x[f]) for x in rows]; z,w=crt(rr,PRIMES),garner(rr,PRIMES); bound=72*fact if f.startswith('E') else fact
   assert z==w and 0<=z<=bound and all(z%p==r for p,r in zip(PRIMES,rr));vals[f[:-8]]=z;glob[f[:-8]]+=m*z
  for t in range(1,6):
   assert vals[f'far3_h{t}']<=vals[f'far2_h{t}']<=vals[f'far1_h{t}']<=vals[f'h{t}']
  support[str(v)]={'far3_le_far2_le_far1_le_h':True};r45[str(v)]=vals['h3']-3*vals['h4']+2*vals['h5'];out[str(v)]={'multiplicity':m,'values':vals};
 for t in range(1,6): assert glob[f'h{t}']==int(old['global']['h'][f'h{t}'])
 A2=int(old['global']['A2']);h={f'h{i}':glob[f'h{i}'] for i in range(1,6)};E={f'E{i}':glob[f'E{i}'] for i in range(1,6)}
 result={'status':'passed','scope':'exact finite n=7 spatial h/E/far CRT reconstruction; no asymptotic conclusion and no per-orbit R45 assertion','source_sha256':SOURCE,'parameters':{'n':7,'orbits':MULT,'multiplicity_sum':49,'primes':PRIMES,'crt_modulus':mod,'h_far_upper_bound':fact,'E_upper_bound':72*fact,'euclidean_squared_maximum':72},'orbits':out,'global':{'A2':A2,'h':h,'E':E,'h_matches_public_n7_exact_all_5':True,'h3_minus_3h4_plus_2h5':h['h3']-3*h['h4']+2*h['h5'],'first_five_sum_t_h_over_A2':frac(sum(t*h[f'h{t}'] for t in range(1,6)),A2),'first_five_sum_t_E_over_49A2':frac(sum(t*E[f'E{t}'] for t in range(1,6)),49*A2),'far':{f'far{k}':{f'h{i}':glob[f'far{k}_h{i}'] for i in range(1,6)} for k in range(1,4)}},'checks':{'Pratt_primes_verified':True,'direct_and_Garner_CRT_agree':True,'all_40_raw_files_valid':True,'D4_multiplicity_sum_49':True,'all_h_match_public_exact':True,'per_orbit_far_support_monotonic':support,'per_orbit_R45_values_recorded_not_assumed_zero':r45}}
 a.output.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n');print(json.dumps({'status':'passed','h':h,'first5_h_over_A2':result['global']['first_five_sum_t_h_over_A2'],'first5_E_over_49A2':result['global']['first_five_sum_t_E_over_49A2']},indent=2))
if __name__=='__main__':main()
