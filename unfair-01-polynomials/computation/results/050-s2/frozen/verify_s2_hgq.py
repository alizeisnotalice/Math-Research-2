"""Independent S2 unrestricted cofactor trees using integer H(GQ) convolution.

H=P*Pbar is integral. Each state stores G,Q and W=GQ; C_k is obtained
without forming A or B. Quadratic floors use exact sign bisection, not the
searcher's square-root floor formula. Every rational row keeps both endpoints.
"""
import argparse,functools,hashlib,json,math,os,time,sys
from collections import Counter
from pathlib import Path

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

@functools.lru_cache(maxsize=65536)
def qa_floor(q,T,D):
    if q==0:return 0
    B=abs(q)*(abs(T)+math.isqrt(D)+2)+1;lo,hi=-B,B
    def nonnegative(mid):
        L=q*T-2*mid
        if L>=0 and q>0:return True
        if L<=0 and q<0:return False
        if L>0:return L*L>q*q*D
        return q*q*D>L*L
    assert nonnegative(lo) and not nonnegative(hi)
    while hi-lo>1:
        mid=(hi+lo)//2
        if nonnegative(mid):lo=mid
        else:hi=mid
    return lo

def choices(p,q,T,D):
    z=p+qa_floor(q,T,D)
    return (-z,1-z) if q==0 else (-z,)

def norm(F,T,N):
    H=[0]*15
    for j in range(15):
        H[j]=sum(F[i]*F[j-i] for i in range(8) if 0<=j-i<8)
        if 0<=j-2<8:H[j]+=T*F[j-2]
    H[4]+=N
    return H

def state_hash(states):
    h=hashlib.sha256()
    for G,Q,W in sorted(states,key=lambda z:(z[0],z[1])):
        h.update(json.dumps([G,Q],separators=(',',':')).encode()+b'\n')
    return h.hexdigest()

def tree(row,depth=96,cap=25000,seconds=8):
    T,N,D,b,c,d,e,f=row
    assert D==T*T-4*N and D>0 and math.isqrt(D)**2!=D
    F=[1,b,0,c,d,e,f,1];H=norm(F,T,N);states=[((1,),(1,),(1,))];trace=[];begin=time.monotonic()
    status='unknown_depth'
    for k in range(1,depth+1):
        nxt=[];reject=Counter(forced_endpoint=0,forced_g3=0,nonbinary_C=0)
        for G,Q,W in states:
            pg=sum(F[j]*G[k-j] for j in range(1,min(k,7)+1));pq=sum(F[j]*Q[k-j] for j in range(1,min(k,7)+1))
            ag=G[k-2] if k>=2 else 0;aq=Q[k-2] if k>=2 else 0;pq+=T*aq
            internal=sum(G[i]*Q[k-i] for i in range(1,k))
            R=internal+sum(H[j]*W[k-j] for j in range(1,min(k,14)+1))
            for g in choices(pg,ag,T,D):
                for q in choices(pq,-aq,T,D):
                    if k==3 and (g!=0 or q!=0):reject['forced_g3']+=1;continue
                    if k==5 and (pg+g!=0 or ag!=0 or pq+q!=0 or aq!=0):reject['forced_endpoint']+=1;continue
                    if R+g+q not in (0,1):reject['nonbinary_C']+=1;continue
                    nxt.append((G+(g,),Q+(q,),W+(internal+g+q,)))
        assert len({(g,q) for g,q,w in nxt})==len(nxt)
        states=nxt;trace.append({'k':k,'live':len(states),'state_sha256':state_hash(states),'rejected':dict(reject)})
        if not states:status='dead_all_prefixes';break
        if len(states)>cap:status='unknown_state_cap';break
        if time.monotonic()-begin>seconds:status='unknown_time_limit';break
    out={'row':row,'status':status,'depth_reached':len(trace),'trace':trace}
    if status!='dead_all_prefixes':out['live_prefixes']=[{'G':list(g),'Q':list(q)} for g,q,w in states]
    return out

def atomic(p,obj):
    tmp=p.with_suffix('.partial')
    with tmp.open('w') as f:json.dump(obj,f,separators=(',',':'));f.write('\n');f.flush();os.fsync(f.fileno())
    os.replace(tmp,p)

def main():
    if not __debug__:raise RuntimeError('no -O')
    ap=argparse.ArgumentParser();ap.add_argument('--screen',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);ap.add_argument('--start',type=int,default=0);ap.add_argument('--count',type=int,default=359);ap.add_argument('--seconds',type=float,default=14400);args=ap.parse_args();args.out.mkdir(parents=True,exist_ok=True)
    assert args.start>=0 and args.start+args.count<=359
    code=sha(Path(__file__));begin=time.monotonic();done=[]
    for i in range(args.start,args.start+args.count):
        if time.monotonic()-begin>=args.seconds:break
        source=args.screen/f'core-{i:03}.json'
        while not source.exists() and time.monotonic()-begin<args.seconds:time.sleep(2)
        if not source.exists():break
        source_hash=sha(source);s=json.loads(source.read_text());assert s['status']=='complete' and s['index']==i
        out=args.out/f'core-{i:03}.json'
        if out.exists():
            x=json.loads(out.read_text());assert x['source_sha256']==code and x['screen_sha256']==source_hash;done.append(i);continue
        records=[tree(row) for row in s['survivors']];counts=Counter(r['status'] for r in records)
        assert len({tuple(r['row']) for r in records})==len(records)
        x={'status':'completed_all_listed_mothers','index':i,'source_sha256':code,'screen_sha256':source_hash,'screen_binding':s['binding'],'mothers':len(records),'status_counts':dict(counts),'records':records}
        atomic(out,x);done.append(i)
        if len(done)%20==0:print(json.dumps({'completed_TN':len(done),'index':i,'elapsed_seconds':time.monotonic()-begin}),flush=True)
    statuses=Counter();mothers=0;max_depth=0;max_live=0
    for i in done:
        x=json.loads((args.out/f'core-{i:03}.json').read_text());statuses.update(x['status_counts']);mothers+=x['mothers']
        for row in x['records']:
            max_depth=max(max_depth,row['depth_reached']);max_live=max(max_live,max(z['live'] for z in row['trace']))
    summary={'status':'complete' if len(done)==args.count else 'partial_unknown','requested':args.count,'completed_TN':len(done),'indices':done,'mothers':mothers,'status_counts':dict(statuses),'max_depth':max_depth,'max_live':max_live,'source_sha256':code,'elapsed_seconds':time.monotonic()-begin,'python':sys.version}
    atomic(args.out/'summary.json',summary);print(json.dumps(summary),flush=True)

if __name__=='__main__':main()
