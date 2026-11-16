"""Independent fixed H10(2,3) exact-certificate checker.

The final certificate is intentionally not assumed to exist at preparation
time.  This checker does not read any options/search/cache file: layouts are
regenerated from 3u+2v=E.  It uses the separately audited, exact helper
modules copied with this review directory.
"""
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
import argparse, hashlib, json, math
from root_ball_common import complex_from_json, complex_from_strings, parse_fraction, polynomial_rouche_bounds, lift_rouche_bounds, disk_disjoint, abs_bounds
from verify_modular_galois import verify as verify_galois

H=[1,1,1,0,0,0,0,0,0,1,1]
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def integer(x): assert type(x) is int; return x
def unique(rows,key):
 out={}
 for r in rows:
  k=key(r);assert k not in out;out[k]=r
 return out
def conv(a,b):return [sum(a[i]*b[k-i] for i in range(len(a)) if 0<=k-i<len(b)) for k in range(len(a)+len(b)-1)]
def stretch(a,n):z=[0]*(n*(len(a)-1)+1);z[::n]=a;return z
def disk(r):
 d=r['disk'];z,q=complex_from_json(d['center']),parse_fraction(d['radius']);assert q>0;return z,q
def separate(ds):
 for i,(z,r) in enumerate(ds):
  assert abs(z[1])>r
  for w,s in ds[:i]:assert disk_disjoint(z,r,w,s)
def coef(factors,ids,k):
 lo=hi=[F(1)]
 for g in sorted(ids):
  n=min(k+1,len(lo)+2);a=[];b=[]
  for t in range(n):
   ls=us=F(0)
   for j in range(max(0,t-len(lo)+1),min(2,t)+1):
    x,y=lo[t-j],hi[t-j];p,q=factors[g][j]
    v=(x*p,x*q,y*p,y*q);ls+=min(v);us+=max(v)
   a.append(ls);b.append(us)
  lo,hi=a,b
 return lo[k],hi[k]
def verify(inp,balls_p,cert_p,gal_p):
 source=json.loads(Path(inp).read_text());assert isinstance(source['candidates'],list) and len(source['candidates'])==1
 raw=source['candidates'][0];balls=json.loads(Path(balls_p).read_text());cert=json.loads(Path(cert_p).read_text());gal=json.loads(Path(gal_p).read_text())
 assert raw['H_coefficients_ascending']==H and (raw['a'],raw['b'],raw['degree'])==(2,3,50) and H!=H[::-1]
 C=conv(stretch(H,2),stretch(H,3));assert set(C)<={0,1} and sum(C)==25 and raw['C_coefficients_ascending']==C
 gr=gal['certificates'];assert len(gr)==1 and gr[0]['coefficients_ascending']==H;gv=verify_galois(gr[0])
 assert balls['H']==H and (balls['a'],balls['b'])==(2,3) and balls['input_sha256']==sha(inp)
 assert cert['format']=='H10-a2b3-negative-Q-compact-v1' and cert['input_sha256']==sha(inp) and cert['balls_sha256']==sha(balls_p)
 assert cert['H']==H and (cert['a'],cert['b'],cert['codes_total'])==(2,3,8310)
 scale=int(balls['scale']);assert scale>0
 groups=unique(raw['groups'],lambda x:integer(x['index']));blocks=unique(raw['blocks'],lambda x:integer(x['block']));moms=unique(balls['mothers'],lambda x:integer(x['block']));lifts=unique(balls['lifts'],lambda x:integer(x['source_group_index']))
 assert set(groups)==set(lifts)==set(range(25)) and set(blocks)==set(moms)==set(range(5))
 md={};ds=[]
 for b,m in moms.items():
  z,r=disk(m);assert m['id']==f'mother-{b}' and z[1]>r and z==complex_from_strings(blocks[b]['parent_root_positive_imag'])
  assert m['kind']=='mother'
  L,R=polynomial_rouche_bounds(H,z,r,scale);assert L<R and parse_fraction(m['lhs'])==L and parse_fraction(m['rhs'])==R
  md[b]=(z,r);ds += [(z,r),((z[0],-z[1]),r)]
 separate(ds);fac={};slot={};td=[]
 for i,g in groups.items():
  b,f,j=integer(g['block']),g['family'],integer(g['lift_index']);q=2 if f=='a_lift' else 3
  assert b in md and f in ('a_lift','b_lift') and 0<=j<q and (b,f,j) not in slot
  x=lifts[i];assert (x['block'],x['family'],x['lift_index'],x['exponent'])==(b,f,j,q) and x['id']==f'{f}-{b}-{j}' and x['parent']==f'mother-{b}' and x['kind']=='lift'
  z,r=disk(x);assert z==complex_from_strings(g['lift_root']) and md[b][0]==complex_from_strings(g['parent_root_positive_imag'])
  L,R=lift_rouche_bounds(*md[b],z,r,scale,q);assert L<R and parse_fraction(x['lhs'])==L and parse_fraction(x['rhs'])==R
  lo,hi=abs_bounds(z,scale);lo-=r;hi+=r;assert lo>0;fac[i]=[(lo*lo,hi*hi),(-2*(z[0]+r),-2*(z[0]-r)),(F(1),F(1))];slot[b,f,j]=i;td += [(z,r),((z[0],-z[1]),r)]
 separate(td);slots=[]
 for b in range(5):
  aa=[slot[b,'a_lift',j] for j in range(2)];bb=[slot[b,'b_lift',j] for j in range(3)]
  assert blocks[b]['a_lift_groups']==aa and blocks[b]['b_lift_groups']==bb and blocks[b]['all_group_indices']==aa+bb;slots.append(aa+bb)
 layouts={};counts={}
 for E in (2,3,4,5,6):
  if E==6:cs=[[(0,1)],[(0,1),(2,3,4)],[(0,1),(2,3,4)],[(0,1),(2,3,4)],[(0,1),(2,3,4)]]
  else:
   uv=[(u,v) for u in range(3) for v in range(4) if 3*u+2*v==E];assert len(uv)==1;u,v=uv[0]
   ch=sorted([tuple(a)+tuple(2+i for i in d) for a in combinations(range(2),u) for d in combinations(range(3),v)],key=lambda t:sum(1<<j for j in t));cs=[ch]*5
  n=math.prod(map(len,cs));counts[E]=n
  for code in range(n):
   z=code;digits=[];picked=[]
   for b,ch in enumerate(cs):d=z%len(ch);z//=len(ch);digits.append(d);picked += [slots[b][j] for j in ch[d]]
   assert z==0 and len(set(picked))==len(picked);layouts[E,code]=(digits,set(picked))
 assert counts=={2:243,3:32,4:243,5:7776,6:16} and len(layouts)==8310
 rec=unique(cert['negative_certificates'],lambda r:(integer(r['E']),integer(r['code'])))
 unk=unique(cert['unknown'],lambda r:(integer(r['E']),integer(r['code'])))
 assert set(unk)=={(6,0)} and len(rec)==8309 and not set(rec)&set(unk) and set(rec)|set(unk)==set(layouts)
 assert unk[6,0]['digits']==[0]*5 and unk[6,0]['terminal_standard'] is True
 sel=layouts[6,0][1];assert sel=={i for i,g in groups.items() if g['family']=='a_lift'} and set(groups)-sel=={i for i,g in groups.items() if g['family']=='b_lift'} and len(sel)*2==20 and (25-len(sel))*2==30
 bits=integer(cert['denominator_power']);assert 1<=bits<=1024;den=1<<bits;margins=[]
 for key,r in rec.items():
  digits,picked=layouts[key];assert r['digits']==digits and r['side'] in ('A','B');ids=picked if r['side']=='A' else set(groups)-picked;k=r['coefficient_index'];assert type(k) is int and 1<=k<2*len(ids)
  l,u=F(int(r['lo_n']),den),F(int(r['hi_n']),den);assert l<=u<0;L,U=coef(fac,ids,k);assert l<=L<=U<=u;margins.append(-u)
 return {'status':'all_8310_layouts_verified_exactly','root_disks':50,'negative_layouts':8309,'standard_layouts':1,'layout_counts':counts,'minimum_certified_negative_margin':str(min(margins)),'Galois':gv,'input_sha256':sha(inp),'balls_sha256':sha(balls_p),'certificate_sha256':sha(cert_p),'galois_sha256':sha(gal_p),'verifier_sha256':sha(__file__),'arithmetic':'integers and Fraction only; no endpoint pinning; no options cache'}
def main():
 p=argparse.ArgumentParser();[p.add_argument('--'+x,type=Path,required=True) for x in ('input','balls','certificate','galois','record')];a=p.parse_args();o=verify(a.input,a.balls,a.certificate,a.galois);a.record.write_text(json.dumps(o,indent=2)+'\n');print(json.dumps(o))
if __name__=='__main__':main()
