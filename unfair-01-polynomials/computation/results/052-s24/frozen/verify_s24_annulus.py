"""Independent complete root-screen replay: trace/norm formula, descending Schur.

Rebuilds deduplicated mother parameters from all row7 states, and compares
every first-failure byte, not only survivor counts. No parent algorithm import.
"""
import argparse,concurrent.futures as cf,gzip,hashlib,json,math,os,time
from pathlib import Path
from collections import Counter
ROOT=Path(__file__).resolve().parent
POINTS=((1,1),(1,2),(2,1),(1,4),(4,1),(3,4),(4,3))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def atom(p,obj):
    q=p.with_suffix('.partial')
    with q.open('w') as f:json.dump(obj,f,separators=(',',':'));f.write('\n');f.flush();os.fsync(f.fileno())
    q.replace(p)
def schur(d):
    while len(d)>1:
        g=math.gcd(*d)
        if d[0]<0:g=-g
        d=[x//g for x in d]
        if abs(d[-1])>=d[0]:return len(d)-1
        first,last=d[0],d[-1];n=len(d)-1;d=[first*d[i]-last*d[n-i] for i in range(n)]
    return None
def coordinates(m,n,D):return (2*m+n,n,2) if D%4==1 else(m,n,1)
def norm(F,ta,td,na,nd,cross):
    H=[0]*15
    for i,u in enumerate(F):
        if not u:continue
        for j,v in enumerate(F):
            if v:H[i+j]+=u*v
        H[i+2]+=ta*u;H[i+4]+=td*u
    H[4]+=na;H[6]+=cross;H[8]+=nd;return H
def one(D,input_dir,parent_dir,out_dir,deadline,source_hash):
    inp=Path(input_dir)/f'field-D{D}.json';par=Path(parent_dir)/f'field-D{D}.json';out=Path(out_dir)/f'field-D{D}.json'
    while not par.exists() and time.time()<deadline:time.sleep(2)
    if not par.exists():return {'D':D,'status':'partial_missing_parent'}
    binding={'source_sha256':source_hash,'input_sha256':sha(inp),'parent_sha256':sha(par)}
    if out.exists():z=json.loads(out.read_text());assert z['binding']==binding and z['status']=='passed';return z
    started=time.monotonic();x=json.loads(inp.read_text());y=json.loads(par.read_text());assert x['D']==y['D']==D and x['status']==y['status']=='complete';mothers={}
    for r in x['rows']:
        key=tuple(r[:4]+r[10:12]);fs=tuple(r[16]);assert fs==tuple(sorted(set(fs)))
        if key in mothers:assert mothers[key]==fs
        else:mothers[key]=fs
    del x
    packed=(Path(parent_dir)/y['ledger_file']).read_bytes();assert hashlib.sha256(packed).hexdigest()==y['ledger_sha256'];codes=gzip.decompress(packed);index=0;sur=[];counts=Counter(unique_bases=len(mothers))
    for (am,an,dm,dn,c,e),fs in sorted(mothers.items()):
        R,V,s=coordinates(am,an,D);U,W,t=coordinates(dm,dn,D);assert s==t
        ta=2*R//s;td=2*U//s;na=(R*R-D*V*V)//(s*s);nd=(U*U-D*W*W)//(s*s);cross=2*(R*U-D*V*W)//(s*s)
        assert s*ta==2*R and s*td==2*U and s*s*na==R*R-D*V*V and s*s*nd==U*U-D*W*W and s*s*cross==2*(R*U-D*V*W)
        tests=[]
        for p,q in POINTS:
            powers=[p**i*q**(7-i) for i in range(8)];real=s*(powers[0]+c*powers[3]+e*powers[5]+powers[7])+R*powers[2]+U*powers[4];imag=V*powers[2]+W*powers[4];tests.append((real,s*powers[6],imag*imag*D))
        for f in fs:
            counts['unique_mothers']+=1;code=None
            for j,(u,v,z) in enumerate(tests):
                real=u+v*f
                if real<=0 or real*real<=z:code=10+j;counts[f'axis{j}']+=1;break
            if code is None:
                counts['positive_axis']+=1;H=norm([1,0,0,c,0,e,f,1],ta,td,na,nd,cross);n=schur([(v<<i) for i,v in enumerate(H)][::-1])
                if n is not None:code=20+n;counts['outer_rejected']+=1
                else:
                    n=schur([H[j]<<(14-j) for j in range(15)])
                    if n is not None:code=40+n;counts['inner_rejected']+=1
                    else:code=255;counts['survived']+=1;sur.append([D,am,an,dm,dn,c,e,f])
            assert index<len(codes) and codes[index]==code,('failure_byte_mismatch',D,index);index+=1
    assert index==len(codes) and sur==y['survivors'] and Counter(y['counts'])==counts
    z={'status':'passed','D':D,'binding':binding,'verified_bytes':index,'raw_codes_sha256':hashlib.sha256(codes).hexdigest(),'counts':dict(counts),'seconds':time.monotonic()-started};atom(out,z);return z
def main():
    if not __debug__:raise RuntimeError('no -O')
    a=argparse.ArgumentParser();a.add_argument('--input',type=Path,default=ROOT/'row7-full');a.add_argument('--parent',type=Path,default=ROOT/'annulus-full');a.add_argument('--out',type=Path,required=True);a.add_argument('--workers',type=int,default=2);a.add_argument('--seconds',type=int,default=14400);p=a.parse_args();p.out.mkdir(parents=True,exist_ok=True);ds=json.loads((p.input/'summary.json').read_text())['fields'];deadline=time.time()+p.seconds;code=sha(Path(__file__));done=[];counts=Counter()
    with cf.ProcessPoolExecutor(max_workers=p.workers) as pool:
        futures=[pool.submit(one,D,str(p.input),str(p.parent),str(p.out),deadline,code) for D in ds]
        for f in cf.as_completed(futures):
            z=f.result();done.append(z);counts.update(z.get('counts',{}));print(json.dumps({k:v for k,v in z.items() if k!='binding'}),flush=True);atom(p.out/'progress.json',{'completed':len(done),'expected':len(ds),'fields':[z['D'] for z in done]})
    out={'status':'passed_all_fields' if all(z['status']=='passed' for z in done) else 'partial_unknown','fields':len(done),'counts':dict(counts),'source_sha256':code,'scope':'Exact necessary root-screen replay only; no cofactor exclusion.'};atom(p.out/'summary.json',out);print(json.dumps(out),flush=True)
if __name__=='__main__':main()
