#!/usr/bin/env python3
"""Exact rational Rouché disks for H8(x^2)H8(x^3), using input suggestions only."""
from fractions import Fraction
from pathlib import Path
import argparse, hashlib, json
from root_ball_common import complex_from_strings, complex_to_json, disk_disjoint, disk_nonreal, fraction_string, lift_rouche_bounds, parse_fraction, polynomial_rouche_bounds

H=[1,1,1,0,0,0,0,0,0,1,1]; SCALE=10**160; MR=Fraction(1,10**70); LR=Fraction(1,10**65)
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def atom(p,x):
 p=Path(p); t=p.with_suffix(p.suffix+'.tmp'); t.write_text(json.dumps(x,indent=2,sort_keys=True)+'\n'); t.replace(p)
def dj(c,r): return {'center':complex_to_json(c),'radius':fraction_string(r)}
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--input',type=Path,required=True);ap.add_argument('--output',type=Path,required=True); a=ap.parse_args()
 raw=json.loads(a.input.read_text())['candidates'][0]
 if raw.get('H_coefficients_ascending')!=H or (raw.get('a'),raw.get('b'),raw.get('degree'))!=(2,3,50): raise SystemExit('not fixed H10 a2b3 input')
 gs=raw['groups']; byblock={}
 for g in gs: byblock.setdefault(int(g['block']),g)
 if set(byblock)!=set(range(5)) or len(gs)!=25: raise SystemExit('expected 5 blocks /25 lift groups')
 moms={}; entries=[]
 for b,g in sorted(byblock.items()):
  c=complex_from_strings(g['parent_root_positive_imag']); L,R=polynomial_rouche_bounds(H,c,MR,SCALE)
  if not L<R or not disk_nonreal(c,MR): raise SystemExit(f'mother {b} fails')
  moms[b]=(c,MR); entries.append({'id':f'mother-{b}','kind':'mother','block':b,'disk':dj(c,MR),'lhs':fraction_string(L),'rhs':fraction_string(R)})
 lifts=[]
 for g in gs:
  b=int(g['block']); fam=str(g['family']); q=2 if fam=='a_lift' else 3
  c=complex_from_strings(g['lift_root']); ac,ar=moms[b];L,R=lift_rouche_bounds(ac,ar,c,LR,SCALE,q)
  if not L<R or not disk_nonreal(c,LR): raise SystemExit(f'lift {g["index"]} fails')
  lifts.append({'id':f'{fam}-{b}-{int(g["lift_index"])}','kind':'lift','family':fam,'exponent':q,'block':b,'lift_index':int(g['lift_index']),'source_group_index':int(g['index']),'parent':f'mother-{b}','disk':dj(c,LR),'lhs':fraction_string(L),'rhs':fraction_string(R)})
 disks=[]
 for x in lifts:
  c=x['disk']['center'];z=(parse_fraction(c['re']),parse_fraction(c['im']));r=parse_fraction(x['disk']['radius']);disks += [(x['id'],z,r),(x['id']+'-conj',(z[0],-z[1]),r)]
 for i,(_,c,r) in enumerate(disks):
  if not disk_nonreal(c,r) or any(not disk_disjoint(c,r,d,s) for _,d,s in disks[:i]): raise SystemExit('lift disks not nonreal/disjoint')
 out={'format':'H10-a2b3-rational-rouche-balls-v1','input_sha256':sha(a.input),'H':H,'a':2,'b':3,'scale':str(SCALE),'mother_radius':fraction_string(MR),'lift_radius':fraction_string(LR),'mothers':entries,'lifts':lifts,'coverage':{'target_conjugate_groups':len(lifts),'target_disks':2*len(lifts),'pairwise_disjoint':True}}
 a.output.parent.mkdir(parents=True,exist_ok=True);atom(a.output,out);print(json.dumps({'status':'generated','groups':len(lifts),'sha256':sha(a.output)}))
if __name__=='__main__': main()
