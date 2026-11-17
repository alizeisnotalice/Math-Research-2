"""Canonical S234 exact initial box and fourth-row necessary prefix screen.

Conjugate exchange permits a_n>0. This halves the earlier two-orientation
box. Integer k scanning never relies on a floating-point range calculation.
Records are compact tuples; no e/f or arbitrary cofactor claim is made.
"""
import argparse,functools,hashlib,importlib.util,json,math,os,sys,time
from pathlib import Path
from collections import Counter
ROOT=Path(__file__).resolve().parent
RING=ROOT/'old_ring.py'
s=importlib.util.spec_from_file_location('ring',RING);m=importlib.util.module_from_spec(s);s.loader.exec_module(m)
def floorq(x,D):
    R,V,den=m.rv(x,D);B=abs(R)+abs(V)*(D+1)+2;lo,hi=-B,B
    while hi-lo>1:
        mid=(lo+hi)//2
        if m.sign(m.add(x,(-mid,0)),D)>=0:lo=mid
        else:hi=mid
    return lo
def unit(x,D):return m.sign(x,D)>=0 and m.sign(m.add((1,0),m.neg(x)),D)>=0
def shifts(x,D):
    n=-floorq(x,D)
    return [j for j in (n-1,n,n+1) if unit(m.add(x,(j,0)),D)]
def lattice(D,num,den,positive):
    M=(2*num+den-1)//den if D%4==1 else (num+den-1)//den
    return [(u,v) for v in range(1 if positive else -M,M+1) if v for u in range(-M,M+1) if m.abs_both_lt((u,v),D,num,den)]
@functools.lru_cache(maxsize=131072)
def row4(D,am,an,dm,dn):
    a=(am,an);d=(dm,dn);T=2*am+an if D%4==1 else 2*am;u2=-floorq(a,D);v2=1-T-u2;A2=m.add(a,(u2,0));B2=m.add(m.conj(a,D),(v2,0))
    U=m.add(d,m.scale(u2,a));V=m.add(m.conj(d,D),m.scale(v2,m.conj(a,D)));base=m.add(m.add(U,V),m.mul(A2,B2,D));assert base[1]==0
    out=[]
    for u4 in shifts(U,D):
        v4=1-base[0]-u4
        if unit(m.add(V,(v4,0)),D):out.append((u2,v2,u4,v4))
    assert len(out)<=1;return out[0] if out else None
def one(D):
    aa=lattice(D,117,4,True);dd=lattice(D,419,8,False);counts=Counter(raw_a=len(aa),raw_d=len(dd));rows=[]
    for ell in [*range(-8,0),*range(1,9)]:
        for am,an in aa:
            a=(am,an)
            if not m.abs_both_lt(m.add((ell*ell,0),m.scale(-2,a)),D,55,4):counts['s2_reject']+=1;continue
            counts['s2_pass']+=1;pa=m.scale(2,m.mul(a,a,D))
            for k in range(-180,181):
                if abs(-ell**3-3*k)*8>=203:continue
                b=m.add(m.scale(ell,a),(k,0))
                if not m.abs_both_lt(b,D,419,8):counts['b_reject']+=1;continue
                counts['ak_pass']+=1;U=pa[0]+ell**4+4*ell*k;V=pa[1]
                for dm,dn in dd:
                    if not m.abs_both_lt((U-4*dm,V-4*dn),D,787,16):counts['s4_reject']+=1;continue
                    counts['core_pass']+=1;low=row4(D,am,an,dm,dn)
                    if low is None:counts['row4_reject']+=1;continue
                    counts['row4_pass']+=1;rows.append([ell,am,an,k,dm,dn,*low])
    return dict(counts),rows
def atomic(p,x):
    tmp=p.with_suffix('.partial')
    with tmp.open('w') as f:json.dump(x,f,separators=(',',':'));f.write('\n');f.flush();os.fsync(f.fileno())
    tmp.replace(p)
def main():
    if not __debug__:raise RuntimeError('no -O')
    ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True);ap.add_argument('--seconds',type=float,default=14400);args=ap.parse_args();args.out.mkdir(parents=True,exist_ok=True)
    fields=[D for D in range(2,190) if (D<=47 or D%4==1) and all(D%(p*p) for p in range(2,math.isqrt(D)+1))];assert len(fields)==58
    binding={'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'ring_sha256':hashlib.sha256(RING.read_bytes()).hexdigest()};start=time.monotonic();done=[];total=Counter()
    for D in fields:
        p=args.out/f'field-D{D}.json'
        if p.exists():x=json.loads(p.read_text());assert x['status']=='complete' and x['binding']==binding and x['D']==D
        else:
            if time.monotonic()-start>args.seconds:break
            counts,rows=one(D);x={'status':'complete','D':D,'binding':binding,'counts':counts,'rows':rows};atomic(p,x)
        done.append(D);total.update(x['counts']);print(json.dumps({'D':D,'completed':len(done),'core_pass':total['core_pass'],'row4_pass':total['row4_pass'],'seconds':time.monotonic()-start}),flush=True)
    if len(done)==58:assert total['ak_pass']==61596 and total['core_pass']==5181392
    out={'status':'complete' if len(done)==58 else 'partial_unknown','fields':done,'counts':dict(total),'binding':binding,'seconds':time.monotonic()-start,'python':sys.version,'record_schema':['ell','a_m','a_n','k','d_m','d_n','u2=G2-ell^2','v2=Q2-ell^2','u4=G4-ell^2*G2-2ell*k','v4=Q4-ell^2*Q2-2ell*k'],'scope':'Canonical a_n>0, complete initial s2/s3/s4 and coefficient box plus row4 only. No e/f or full cofactor exclusion.'};atomic(args.out/'summary.json',out);print(json.dumps(out),flush=True)
if __name__=='__main__':main()
