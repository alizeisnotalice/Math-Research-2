"""Independent integer convolution checks; not a proof by finite sampling."""
import argparse, hashlib, json, math, random, sys, time
from pathlib import Path

def mul(a,b):
    c=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b): c[i+j]+=x*y
    return c
def stretch(a,k):
    b=[0]*((len(a)-1)*k+1)
    for i,x in enumerate(a): b[k*i]=x
    return b
def at(a,i):return a[i] if 0<=i<len(a) else 0
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--record',type=Path,required=True);a=ap.parse_args()
    assert Path('/Volumes/data').is_mount()
    start=time.monotonic();rng=random.Random(2026111633);count=0;pairs=[]
    for p in range(2,10):
        for q in range(p+1,15):
            if math.gcd(p,q)!=1:continue
            pairs.append([p,q]);F=[rng.randrange(-4,5) for _ in range(6*q+1)]
            G=[rng.randrange(0,5) for _ in range(2*q)];F[0]=G[0]=1
            A=mul(stretch(F,p),stretch(G,q));B=mul(stretch(F,q),stretch(G,p))
            for r in range(len(F)):
                assert at(A,p*r)==sum(at(F,r-q*k)*at(G,p*k) for k in range(r//q+1))
                assert at(B,q*r)==at(F,r)+at(F,r-p)*at(G,q)
                count+=2
                for s in range(1,q):
                    assert at(B,q*r+p*s)==at(F,r)*at(G,s)+at(F,r-p)*at(G,s+q)
                    assert at(B,q*r+p*(s+q))==at(F,r)*at(G,s+q)+at(F,r+p)*at(G,s)
                    count+=2
    branches=[]
    for p,q in pairs:
        if q<2*p or (p,q)==(3,5):
            cap=2*q if (p,q)==(3,5) else min(2*q,3*p)
            for ell in range(2,(cap-1)//p+1):
                r=q*ell
                if ell==2:
                    s=2*p-q;assert 0<s<q
                    v=r-q;w=r+p-q
                    assert (r-2*q)+q==v and (r-p)+s==w
                else:
                    assert (p,q,ell)==(3,5,3)
                    v=r-6;w=r-3
                    assert (r-15)+9==v
                assert 0<=v<w<r and w-v==p
                assert q*v+p*q==q*w
                H=[0]*(max(w,q)+1)
                for k in [0,q,v,w]:H[k]=1
                C=mul(stretch(H,p),stretch(H,q));assert C[q*w]>=2
                branches.append({'p':p,'q':q,'ell':ell,'r':r,'collision_index':q*w,'value':C[q*w]})
    rec={'status':'passed','python':sys.version,'seed':2026111633,'pairs':pairs,'coefficient_identity_checks':count,'branches':branches,'elapsed_seconds':time.monotonic()-start,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'scope':'Finite exact integer formula checks only; general theorem requires the analytic proof.'}
    a.record.write_text(json.dumps(rec,indent=2)+'\n');print(json.dumps({k:rec[k] for k in ['status','coefficient_identity_checks','elapsed_seconds']}))
if __name__=='__main__':main()
