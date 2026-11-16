"""Exact small checks supporting, but not replacing, the analytic proof."""
import argparse,hashlib,json,math,random,sys,time
from pathlib import Path
def at(a,k):return a[k] if 0<=k<len(a) else 0
def mul(a,b):
    c=[0]*(len(a)+len(b)-1)
    for i,u in enumerate(a):
        for j,v in enumerate(b):c[i+j]+=u*v
    return c
def poly(indices):
    z=[0]*(max(indices)+1)
    for i,c in indices.items():z[i]=c
    return z
def stretch(a,p):
    z=[0]*((len(a)-1)*p+1)
    for i,c in enumerate(a):z[p*i]=c
    return z
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--record',type=Path,required=True);a=ap.parse_args();start=time.monotonic();rng=random.Random(2026111634)
    isolated=[]
    for p in range(2,101):
        for q in range(p+1,2*p):
            if math.gcd(p,q)!=1:continue
            possible=[j for j in range(2*q) if (j-q)%p==0 and j!=q]
            assert set(possible)<=set([q-p,q+p]);isolated.append([p,q])
    fixed=0
    for Gsupport in [(0,1,5,6),(0,4,5,9)]:
        G=poly({k:1 for k in Gsupport});F=[rng.randrange(-3,4) for _ in range(51)]
        A=mul(stretch(F,2),stretch(G,5));B=mul(stretch(F,5),stretch(G,2))
        shift=15 if 1 in Gsupport else 10
        for i in range(80):
            assert at(A,2*i)==at(F,i)+at(F,i-shift)
            assert at(B,5*i)==at(F,i)+at(F,i-2)
            fixed+=2
    reverse=0
    for q in range(3,11):
        for mask in range(1<<(q-1)):
            T={0}|{i+1 for i in range(q-1) if mask>>i&1};t=max(T);S=T|{q+i for i in T}
            assert {q+t-i for i in S}=={t-i for i in T}|{q+t-i for i in T};reverse+=1
    F=poly({0:1,8:1,9:1,12:-1,15:1,16:1});G=poly({0:1,4:1})
    H=mul(F,G);Y=mul(F,poly({0:1,3:1}));HS={i for i,c in enumerate(H) if c}
    assert H==poly({i:1 for i in [0,4,8,9,13,15,19,20]})
    assert Y==poly({i:1 for i in [0,3,8,9,11,16,18,19]})
    assert not HS&{i+3 for i in HS} and F[12]==-1
    out={'status':'passed','python':sys.version,'seed':2026111634,'singleton_parameter_pairs':len(isolated),'fixed_G_coefficient_checks':fixed,'reverse_support_checks':reverse,'abstract_counterexample_exact':True,'elapsed_seconds':time.monotonic()-start,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'scope':'Finite exact formula checks. The general conclusions rely on the analytic proof.'}
    a.record.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out))
if __name__=='__main__':main()
