"""Complete necessary screens for degree-seven S2; integer arithmetic only.

Input is the independently checked complete low-prefix ledger. Each distinct
(T,N,D,b,c) is expanded over the FULL d/e/f coefficient box. A byte for every
raw mother records its first necessary failure, or 255 for a survivor.
Survival is not a factorization. Arbitrary cofactor trees are separate.
"""
import argparse,gzip,hashlib,json,math,os,sys,time
from collections import Counter
from pathlib import Path
POINTS=((1,1),(1,2),(2,1),(1,4),(4,1),(3,4),(4,3))

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def bound(p,q,T,D,num,den):
    R=2*num-den*abs(2*p+q*T)
    return R>0 and den*den*q*q*D<R*R
def cores(prefix):
    x=json.loads(prefix.read_text());assert len(x['core_records'])==4913 and len(x['survivors'])==395
    out=sorted({tuple(z[k] for k in ('T','N','Delta','b','c')) for z in x['survivors']})
    assert len(out)==359;return out
def disk_failure(a):
    a=list(a)
    while len(a)>1:
        div=math.gcd(*a)
        if a[-1]<0:div=-div
        a=[v//div for v in a]
        if abs(a[0])>=a[-1]:return len(a)-1
        n=len(a)-1;first,last=a[0],a[-1]
        a=[last*a[i]-first*a[n-i] for i in range(1,n+1)]
    return None
def norm(F,T,N):
    H=[0]*15
    for i,v in enumerate(F):
        for j,w in enumerate(F):H[i+j]+=v*w
        H[i+2]+=T*v
    H[4]+=N;return H
def powers(T,N,b,c,d,e,f):
    return (
      (-2*N+b**4+4*b*c-4*d,2*T-4*b*b),
      (5*N*b-b**5-5*b*b*c+5*b*d-5*e,-5*T*b+5*b**3+5*c),
      (2*N*T-9*N*b*b+b**6+6*b**3*c-6*b*b*d+6*b*e+3*c*c-6*f,2*N-2*T*T+9*T*b*b-6*b**4-12*b*c+6*d),
      (-7*N*T*b+14*N*b**3+7*N*c-b**7-7*b**4*c+7*b**3*d-7*b*b*e-7*b*c*c+7*b*f+7*c*d-7,-7*N*b+7*T*T*b-14*T*b**3-7*T*c+7*b**5+21*b*b*c-14*b*d+7*e))
def one(core):
    T,N,D,b,c=core;counts=Counter(raw=105*59*17);ledger=bytearray();sur=[]
    tests=[]
    for p,q in POINTS:
        w=[p**i*q**(7-i) for i in range(8)]
        tests.append((2*(w[0]+b*w[1]+c*w[3]+w[7])+T*w[2],2*w[4],2*w[5],2*w[6],D*w[2]**2))
    for d in range(-52,53):
        p4,q4=powers(T,N,b,c,d,0,0)[0]
        if not bound(p4,q4,T,D,787,16):ledger.extend(bytes([4])*(59*17));continue
        counts['pass4']+=59*17
        for e in range(-29,30):
            pws=powers(T,N,b,c,d,e,0);p5,q5=pws[1]
            if not bound(p5,q5,T,D,3107,32):ledger.extend(bytes([5])*17);continue
            counts['pass5']+=17;p6,q6=pws[2];p7,q7=pws[3]
            for f in range(-8,9):
                if not bound(p6-6*f,q6,T,D,12355,64):ledger.append(6);continue
                counts['pass6']+=1
                if not bound(p7+7*b*f,q7,T,D,49283,128):ledger.append(7);continue
                counts['pass7']+=1;code=None
                for j,(base,dd,ee,ff,V) in enumerate(tests):
                    L=base+dd*d+ee*e+ff*f
                    if L<=0 or L*L<=V:code=10+j;counts[f'axis{j}']+=1;break
                if code is None:
                    counts['positive_axis']+=1;H=norm([1,b,0,c,d,e,f,1],T,N)
                    deg=disk_failure([v<<i for i,v in enumerate(H)])
                    if deg is not None:code=20+deg;counts['outer_rejected']+=1
                    else:
                        deg=disk_failure([v<<i for i,v in enumerate(H[::-1])])
                        if deg is not None:code=40+deg;counts['inner_rejected']+=1
                        else:code=255;counts['annulus_survived']+=1;sur.append([T,N,D,b,c,d,e,f])
                ledger.append(code)
    assert len(ledger)==counts['raw'];return dict(counts),bytes(ledger),sur
def atomic(p,data):
    tmp=p.with_suffix(p.suffix+'.partial')
    with tmp.open('wb') as f:f.write(data);f.flush();os.fsync(f.fileno())
    os.replace(tmp,p)
def main():
    if not __debug__:raise RuntimeError('no -O')
    ap=argparse.ArgumentParser();ap.add_argument('--prefix',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);ap.add_argument('--start',type=int,default=0);ap.add_argument('--count',type=int,default=359);ap.add_argument('--seconds',type=float,default=14400);args=ap.parse_args()
    args.out.mkdir(exist_ok=True,parents=True);allcores=cores(args.prefix);assert 0<=args.start<359 and args.start+args.count<=359
    binding={'source':sha(Path(__file__)),'prefix':sha(args.prefix)};begin=time.monotonic();done=[]
    for i in range(args.start,args.start+args.count):
        p=args.out/f'core-{i:03}.json';z=args.out/f'core-{i:03}-codes.bin.gz'
        if p.exists():
            x=json.loads(p.read_text());assert x['binding']==binding and x['status']=='complete' and x['core']==list(allcores[i]) and x['ledger_sha256']==sha(z);done.append(i);continue
        if time.monotonic()-begin>args.seconds:break
        now=time.monotonic();counts,codes,sur=one(allcores[i]);packed=gzip.compress(codes,mtime=0);atomic(z,packed)
        x={'status':'complete','index':i,'core':allcores[i],'counts':counts,'survivors':sur,'binding':binding,'ledger_sha256':hashlib.sha256(packed).hexdigest(),'raw_code_sha256':hashlib.sha256(codes).hexdigest(),'elapsed_seconds':time.monotonic()-now}
        atomic(p,(json.dumps(x,sort_keys=True)+'\n').encode());done.append(i)
        if len(done)%20==0:print(json.dumps({'completed':len(done),'elapsed_seconds':time.monotonic()-begin}),flush=True)
    counts=Counter();sur=[]
    for i in done:
        x=json.loads((args.out/f'core-{i:03}.json').read_text());counts.update(x['counts']);sur.extend(x['survivors'])
    result={'status':'complete' if len(done)==args.count else 'partial_unknown','requested':args.count,'completed':len(done),'indices':done,'binding':binding,'counts':dict(counts),'survivors':sur,'elapsed_seconds':time.monotonic()-begin,'python':sys.version,'code_schema':'All raw d/e/f lexicographic.4..7=first Newton;10..16=axis;20+degree=outer;40+degree=inner;255=survives. Necessary tests only.'}
    atomic(args.out/'summary.json',(json.dumps(result,sort_keys=True)+'\n').encode());print(json.dumps({k:result[k] for k in ('status','completed','counts','elapsed_seconds')}),flush=True)
if __name__=='__main__':main()
