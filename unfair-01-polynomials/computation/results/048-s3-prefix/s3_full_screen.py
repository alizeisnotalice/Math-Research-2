"""Complete degree-seven S3 necessary screens, checkpointed per trace/norm.

All arithmetic is integer. A surviving mother is NOT a counterexample and is
NOT excluded. Cofactor prefix continuation is a separate task.
"""
import argparse,gzip,hashlib,json,math,os,time,sys
from collections import Counter
from pathlib import Path
POINTS=[(1,1),(1,2),(2,1),(1,4),(4,1),(3,4),(4,3)]

def bound(p,q,T,D,num,den):
    R=2*num-den*abs(2*p+q*T)
    return R>0 and den*den*q*q*D<R*R

def triples():
    out=[]
    for T in range(-16,17):
        for N in range(-71,72):
            D=T*T-4*N
            if D>0 and math.isqrt(D)**2!=D and bound(0,1,T,D,203,24):out.append((T,N,D))
    assert len(out)==669
    return out

def disk_failure(ascending):
    # Schur degree descent on ascending coefficients. GCD normalization is
    # positive up to a common sign chosen to make the leading term positive.
    a=list(ascending)
    while len(a)>1:
        divisor=math.gcd(*a)
        if a[-1]<0:divisor=-divisor
        a=[x//divisor for x in a]
        if abs(a[0])>=a[-1]:return len(a)-1
        first,last=a[0],a[-1];n=len(a)-1
        a=[last*a[i]-first*a[n-i] for i in range(1,n+1)]
    return None

def norm(c,d,e,f,T,N):
    F=[1,0,c,0,d,e,f,1];H=[0]*15
    for i,x in enumerate(F):
        for j,y in enumerate(F):H[i+j]+=x*y
        H[i+3]+=T*x
    H[6]+=N
    return H

def one(triple,allow):
    T,N,D=triple;count=Counter(TN=1);ledger=bytearray();survivors=[]
    tests=[]
    for p,q in POINTS:
        w=[p**i*q**(7-i) for i in range(8)]
        tests.append((2*(w[0]+w[7])+T*w[3],2*w[2],2*w[4],2*w[5],2*w[6],D*w[3]*w[3]))
    for c in range(-6,7):
        count['c']+=1
        es=[e for e in range(-29,30) if bound(-5*e,5*c,T,D,3107,32)]
        permitted=(T,N,D,c) in allow
        for d in range(-52,53):
            if abs(2*c*c-4*d)*16>=787:continue
            count['d']+=1;count['e']+=len(es)
            good_e=[e for e in es if bound(7*(c*e-1),7*(d-c*c),T,D,49283,128)]
            count['s7']+=len(good_e)
            fs=[f for f in range(-8,9) if bound(-3*N-2*c**3+6*c*d-6*f,3*T,T,D,12355,64)]
            count['f']+=len(good_e)*len(fs)
            if not permitted:
                ledger.extend(bytes([0])*(len(good_e)*len(fs)));count['prefix_rejected']+=len(good_e)*len(fs);continue
            for e in good_e:
                for f in fs:
                    count['prefix_allowed']+=1;code=None
                    for j,(L0,Lc,Ld,Le,Lf,V) in enumerate(tests):
                        L=L0+Lc*c+Ld*d+Le*e+Lf*f
                        if L<=0 or L*L<=V:code=1+j;count[f'axis_{j}']+=1;break
                    if code is None:
                        count['positive_axis']+=1;H=norm(c,d,e,f,T,N)
                        failure=disk_failure([v<<i for i,v in enumerate(H)])
                        if failure is not None:code=20+failure;count['outer_rejected']+=1
                        else:
                            failure=disk_failure([v<<i for i,v in enumerate(H[::-1])])
                            if failure is not None:code=40+failure;count['inner_rejected']+=1
                            else:code=255;count['annulus_survived']+=1;survivors.append([T,N,D,c,d,e,f])
                    ledger.append(code)
    assert len(ledger)==count['f']
    return dict(count),bytes(ledger),survivors

def atomic(path,b):
    temp=path.with_suffix(path.suffix+'.partial')
    with temp.open('wb') as f:f.write(b);f.flush();os.fsync(f.fileno())
    os.replace(temp,path)

def main():
    if not __debug__:raise RuntimeError('no -O')
    p=argparse.ArgumentParser();p.add_argument('--prefix',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--start',type=int,default=0);p.add_argument('--count',type=int,default=669);p.add_argument('--seconds',type=float,default=14400);a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True)
    sha=lambda b:hashlib.sha256(b).hexdigest();prefix=json.loads(a.prefix.read_text());assert len(prefix['states'])==17760
    allow={tuple(x[:4]) for x in prefix['states']};binding={'source':sha(Path(__file__).read_bytes()),'prefix':sha(a.prefix.read_bytes())};start=time.monotonic();done=[]
    all_tn=triples();assert a.start>=0 and a.start+a.count<=len(all_tn)
    for i in range(a.start,a.start+a.count):
        target=a.out/f'tn-{i:03}.json';ledger=a.out/f'tn-{i:03}-codes.bin.gz'
        if target.exists():
            x=json.loads(target.read_text());assert x['binding']==binding and x['triple']==list(all_tn[i]) and x['status']=='complete' and x['ledger_sha256']==sha(ledger.read_bytes());done.append(i);continue
        if time.monotonic()-start>a.seconds:break
        begin=time.monotonic();counts,codes,survivors=one(all_tn[i],allow);packed=gzip.compress(codes,mtime=0);atomic(ledger,packed)
        x={'status':'complete','index':i,'triple':all_tn[i],'binding':binding,'counts':counts,'survivors':survivors,'ledger_sha256':sha(packed),'raw_code_sha256':sha(codes),'elapsed_seconds':time.monotonic()-begin}
        atomic(target,(json.dumps(x,sort_keys=True)+'\n').encode());done.append(i)
        if len(done)%10==0:print(json.dumps({'completed':len(done),'index':i,'elapsed':time.monotonic()-start}),flush=True)
    aggregate=Counter();survivors=[]
    for i in done:
        x=json.loads((a.out/f'tn-{i:03}.json').read_text());aggregate.update(x['counts']);survivors.extend(x['survivors'])
    result={'status':'complete' if len(done)==a.count else 'partial_unknown','scope':'Necessary screens only; annulus survivors require cofactor continuation.','start':a.start,'requested':a.count,'completed':len(done),'indices':done,'binding':binding,'counts':dict(aggregate),'survivors':survivors,'elapsed_seconds':time.monotonic()-start,'python':sys.version,'code_schema':'Per T/N, c,d,e,f ascending after Newton bounds. 0=prefix;1..7=first axis;20+degree=outer Schur;40+degree=inner Schur;255=survivor.'}
    atomic(a.out/'summary.json',(json.dumps(result,sort_keys=True)+'\n').encode());print(json.dumps({k:result[k] for k in ['status','completed','counts','elapsed_seconds']}),flush=True)

if __name__=='__main__':main()
