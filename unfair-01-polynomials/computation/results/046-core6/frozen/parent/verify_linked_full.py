#!/usr/bin/env python3
"""Independent full integer replay with one-byte first-failure ledger per input."""
import argparse,hashlib,json,math,time,gzip
from pathlib import Path
HERE=Path(__file__).resolve().parent
SEARCH=HERE.parent/'agents/coprime_stretches/round046/core6_linked'
POINTS=((1,4),(1,2),(3,4),(1,1),(4,3),(2,1),(4,1))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def disk_failure(h):
 a=list(reversed(h))
 while len(a)>1:
  g=math.gcd(*a)
  if a[0]<0:g=-g
  a=[v//g for v in a]
  if abs(a[-1])>=a[0]:return len(a)-1,a[-1],a[0]
  lead,last=a[0],a[-1]
  a=[lead*a[i]-last*a[-i-1] for i in range(len(a)-1)]
 return None
def norm_zero_e(b,t,n,ell,d):
 F=[1,b,0,ell,d,0,1];K=[0,0,1,b,0,0,0];h=[0]*13
 for i in range(7):
  for j in range(7):h[i+j]+=F[i]*F[j]+t*F[i]*K[j]+n*K[i]*K[j]
 return h,[2*F[i]+t*K[i] for i in range(7)]
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--b',type=int,required=True);z=ap.parse_args()
 src=SEARCH/'core6_linked_conditional_ledger_corrected.json'
 all_rows=json.loads(src.read_text())['rows']
 independent={tuple(x) for x in json.loads((HERE/'linked-input-tuples.json').read_text())}
 assert len(all_rows)==len(independent)==35488
 assert {(r['b'],r['T'],r['N'],r['ell']) for r in all_rows}==independent
 rows=[r for r in all_rows if r['b']==z.b];assert rows
 start=time.monotonic();survivors=[];codes=bytearray();digest=hashlib.sha256()
 counts={'mothers':len(rows),'raw_d_e':0,'positive_axis':0,'annulus_schur':0}
 axis_counts={f'{p}/{q}':0 for p,q in POINTS};schur_counts={'H(2x)':0,'reverse(H)(2x)':0}
 for row in rows:
  b,t,n,delta,ell=[row[k] for k in ('b','T','N','delta','ell')]
  tests=[]
  for p,q in POINTS:
   w=[p**i*q**(6-i) for i in range(7)];K=w[2]+b*w[3]
   L=2*(w[0]+b*w[1]+ell*w[3]+w[6])+t*K
   tests.append((L,2*w[4],2*w[5],delta*K*K))
  for d in range(-21,22):
   h0=correction=None
   for e in range(-7,8):
    counts['raw_d_e']+=1;bad=None;code=None
    for index,(L0,Ld,Le,D) in enumerate(tests):
     L=L0+Ld*d+Le*e;margin=L*L-D
     if L<=0 or margin<=0:
      p,q=POINTS[index];bad={'point':[p,q],'L':L,'product_numerator':margin};code=index
      axis_counts[f'{p}/{q}']+=1;kind='axis';break
    if bad is None:
     counts['positive_axis']+=1
     if h0 is None:h0,correction=norm_zero_e(b,t,n,ell,d)
     h=h0.copy()
     for i,c in enumerate(correction):h[i+5]+=e*c
     h[10]+=e*e
     for orientation,coef,offset in [('H(2x)',h,10),('reverse(H)(2x)',h[::-1],30)]:
      ans=disk_failure([v*2**i for i,v in enumerate(coef)])
      if ans is not None:
       deg,constant,leading=ans;bad={'orientation':orientation,'degree':deg,'constant':constant,'leading':leading}
       code=offset+deg;schur_counts[orientation]+=1;kind='schur';break
    if bad is None:
     code=255;counts['annulus_schur']+=1;survivors.append({**row,'d':d,'e':e})
     entry=(b,t,n,delta,ell,d,e,'survivor')
    else:entry=(b,t,n,delta,ell,d,e,kind,bad)
    codes.append(code);digest.update(repr(entry).encode()+b'\n')
  if len(codes)%100000==0:print({'b':b,'checked':len(codes)},flush=True)
 name=str(z.b).replace('-','neg');other=json.loads((SEARCH/f'core6_de_screen_b_{name}.json').read_text())
 assert counts==other['counts'] and survivors==other['survivors']
 assert axis_counts==other['axis_first_failure'] and schur_counts==other['schur_first_failure']
 assert digest.hexdigest()==other['first_failure_stream_sha256']
 ledger=HERE/f'linked-b-{name}-failure-codes.bin.gz';ledger.write_bytes(gzip.compress(bytes(codes),mtime=0))
 out={'status':'passed_complete_independent_b_layer','b':z.b,'counts':counts,'survivors':survivors,'first_failure_stream_sha256':digest.hexdigest(),'compact_code_sha256':sha(ledger),'input_sha256':sha(src),'search_report_sha256':sha(SEARCH/f'core6_de_screen_b_{name}.json'),'verifier_sha256':sha(Path(__file__)),'seconds':time.monotonic()-start,'code_schema':'0..6: first bad positive point; 10+degree: outer Schur;30+degree: reversed Schur;255: necessary survivor. Order is source rows for b, d ascending -21..21, e ascending -7..7.'}
 (HERE/f'linked-b-{name}-summary.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items() if k!='survivors'}),flush=True)
if __name__=='__main__':main()
