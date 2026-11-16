"""Full independent H7 certificate check: exact rational coefficient intervals.
Search code, search group intervals and cached success flags are never imported.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse,json,time,sys
from verify_h7_frontend import verify,H,sha,ALL,require
from verify_modular_galois import verify as modular_verify

def coeff(factors,indices,k):
 lo=[F(1)];hi=[F(1)]
 for i in indices:
  q=factors[i];size=min(k+1,len(lo)+len(q)-1);L=[F(0)]*size;U=[F(0)]*size
  for a in range(len(lo)):
   for j,(c,d) in enumerate(q):
    if a+j>=size:break
    v=(lo[a]*c,lo[a]*d,hi[a]*c,hi[a]*d);L[a+j]+=min(v);U[a+j]+=max(v)
  lo,hi=L,U
 return lo[k],hi[k]
def check(balls,ledger,cert,galois,record):
 b=json.loads(balls.read_text());l=json.loads(ledger.read_text());d,layouts,f=verify(b,l);g=json.loads(galois.read_text());require(len(g['certificates'])==1,'one mother certificate');g=g['certificates'][0];require(g['coefficients_ascending']==H and H!=H[::-1],'S7 fixed nonreciprocal H');gg=modular_verify(g)
 c=json.loads(cert.read_text());require(c['status']=='complete' and c['processed']==c['requested']==12158,'full scan status');require(c['balls_sha256']==sha(balls),'balls binding');bits=c['bits'];require(type(bits)is int and 32<=bits<=512,'valid dyadic bits');D=1<<bits
 standard=sum(1<<x['source_group_index'] for x in b['lifts'] if x['family']=='a');need=set(layouts);seen=set();negative=c['negative_certificates'];unknown=c['unknown'];require(len(negative)==12157 and len(unknown)==1,'record counts')
 for row in negative+unknown:
  E=row['E'];mask=row['mask'];require(type(E)is int and type(mask)is int and (E,mask) in need and (E,mask) not in seen,'unique complete mask records');seen.add((E,mask))
 require(seen==need,'complete layout coverage');u=unknown[0];require(u['E']==12 and u['mask']==standard and u['standard']is True,'exact standard root set')
 require(sum(len(f[i])-1 for i in range(25) if standard>>i&1)==21,'standard a degree21');require(sum(len(f[i])-1 for i in range(25) if (ALL^standard)>>i&1)==28,'standard b degree28')
 start=time.monotonic();minimum=None;Q=1<<96
 for n,row in enumerate(negative,1):
  require(row['side'] in ('A','B'),'side');mask=row['mask'] if row['side']=='A' else ALL^row['mask'];ids=[i for i in range(25) if mask>>i&1];degree=sum(len(f[i])-1 for i in ids);k=row['index'];require(type(k)is int and 0<=k<=degree,'coefficient index')
  a=row['lo_n'];z=row['hi_n'];require(type(a)is str and type(z)is str and len(a)<2000 and len(z)<2000,'integer bound encoding');a=int(a);z=int(z);require(a<=z<0,'proposed negative interval')
  lo,hi=coeff(f,ids,k);require(hi<0,'independent exact negative coefficient')
  # The proposed bounds are suggestions only. Record an independent outward
  # rounded margin from this checker's rational box, never the searcher's bound.
  upper=-((-hi.numerator*Q)//hi.denominator);require(upper<0,'96bit certified margin');margin=-upper;minimum=margin if minimum is None else min(minimum,margin)
  if n%500==0:
   record.with_suffix('.progress.json').write_text(json.dumps({'status':'verifying','checked':n,'elapsed_seconds':time.monotonic()-start})+'\n')
 out={'status':'all_12158_layouts_independently_verified','negative_layouts':len(negative),'standard_layouts':1,'exact_standard_mask':standard,'root_disks':49,'layout_sha256':d['canonical_E_mask_sha256'],'minimum_certified_negative_margin':str(F(minimum,Q)),'S7':gg,'arithmetic':'Exact Fraction box convolution, truncated only above requested coefficient; bound recomputed independently','elapsed_seconds':time.monotonic()-start,'balls_sha256':sha(balls),'ledger_sha256':sha(ledger),'certificate_sha256':sha(cert),'galois_sha256':sha(galois),'verifier_sha256':sha(__file__),'frontend_sha256':sha(Path(__file__).with_name('verify_h7_frontend.py'))}
 record.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out))
if __name__=='__main__':
 sys.set_int_max_str_digits(20000);p=argparse.ArgumentParser()
 for key in ('balls','ledger','certificate','galois','record'):p.add_argument('--'+key,type=Path,required=True)
 a=p.parse_args();check(a.balls,a.ledger,a.certificate,a.galois,a.record)
