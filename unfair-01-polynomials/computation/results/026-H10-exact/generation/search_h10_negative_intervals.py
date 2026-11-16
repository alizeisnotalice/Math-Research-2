#!/usr/bin/env python3
"""Exhaustive rational interval negative-coefficient screen for H8 a2b3 layouts."""
from fractions import Fraction
import sys
from math import prod
from pathlib import Path
import argparse, hashlib, json
from root_ball_common import abs_bounds, complex_from_json, parse_fraction
if hasattr(sys, 'set_int_max_str_digits'): sys.set_int_max_str_digits(100000)

def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def fs(x): return f'{x.numerator}/{x.denominator}'
def atom(p,x):
 p=Path(p);t=p.with_suffix(p.suffix+'.tmp');t.write_text(json.dumps(x,indent=2,sort_keys=True)+'\n');t.replace(p)
def mul(x,y):
 z=[(Fraction(0),Fraction(0)) for _ in range(len(x)+len(y)-1)]
 for i,(a,b) in enumerate(x):
  for j,(c,d) in enumerate(y):
   v=(a*c,a*d,b*c,b*d);lo,hi=z[i+j];z[i+j]=(lo+min(v),hi+max(v))
 return z
def factor(disk,scale):
 c=complex_from_json(disk['center']);r=parse_fraction(disk['radius']);lo,hi=abs_bounds(c,scale)
 cl=max(Fraction(0),lo-r);cu=hi+r
 # ascending (z-root)(z-conjugate): [|root|^2,-2 Re(root),1]
 return [(cl*cl,cu*cu),(-2*(c[0]+r),-2*(c[0]-r)),(Fraction(1),Fraction(1))]
def digits(code,rad):
 a=[]
 for r in rad:a.append(code%r);code//=r
 return a
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--input',type=Path,required=True);ap.add_argument('--options',type=Path,required=True);ap.add_argument('--balls',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args()
 raw=json.loads(a.input.read_text())['candidates'][0];opt=json.loads(a.options.read_text());ball=json.loads(a.balls.read_text())
 if raw['H_coefficients_ascending']!=[1,1,1,0,0,0,0,0,0,1,1] or (opt['a'],opt['b'],opt['block_count'])!=(2,3,5):raise SystemExit('wrong fixed H10 schema')
 if ball['input_sha256']!=sha(a.input):raise SystemExit('ball/input mismatch')
 scale=int(ball['scale']); disks={x['source_group_index']:x['disk'] for x in ball['lifts']}; allidx=set(range(len(raw['groups'])))
 if set(disks)!=allidx:raise SystemExit('ball source coverage mismatch')
 fac={i:factor(disks[i],scale) for i in allidx}; out=[]; unknown=[]; terminal=6
 for Etext in sorted(opt['blocks'][0]['options'],key=int):
  E=int(Etext); choices=[b['options'][Etext] for b in opt['blocks']]
  if E==terminal: choices[0]=[choices[0][0]]
  rad=[len(x) for x in choices]
  for code in range(prod(rad)):
   ds=digits(code,rad); pick=set();
   for ch,d in zip(choices,ds):pick.update(ch[d]['selected_group_indices'])
   for side,ids in (('A',pick),('B',allidx-pick)):
    q=[(Fraction(1),Fraction(1))]
    for i in sorted(ids):q=mul(q,fac[i])
    bad=next(((k,I) for k,I in enumerate(q[1:-1],1) if I[1]<0),None)
    if bad:
     k,I=bad;out.append({'E':E,'code':code,'side':side,'coefficient_index':k,'upper_bound':fs(I[1]),'lower_bound':fs(I[0]),'digits':ds})
     break
   else: unknown.append({'E':E,'code':code,'digits':ds,'terminal_standard':E==terminal and code==0})
 total=len(out)+len(unknown)
 result={'format':'H10-a2b3-negative-Q-interval-v1','input_sha256':sha(a.input),'options_sha256':sha(a.options),'balls_sha256':sha(a.balls),'scale':str(scale),'codes_total':total,'negative_certificates':out,'unknown':unknown,'claim':'Each negative record is for an unnormalised monic selected-root product Q; Q(0)>0, hence Q/Q(0) retains its strict negative coefficient.'}
 atom(a.output,result);print(json.dumps({'negative':len(out),'unknown':len(unknown),'sha':sha(a.output)}))
if __name__=='__main__':main()
