"""S234: Deduplicate row7 prefix mothers and apply exact necessary root screens.

Per field, a 17-bit f mask unions all prefix choices without counting them as
distinct mothers. First-failure byte ledger follows sorted(ell,a,k,d,e), then f.
Four independent field workers; 4h batch. Survivors are not counterexamples.
"""
import argparse,concurrent.futures as cf,gzip,hashlib,importlib.util,json,math,os,sys,time
from pathlib import Path
from collections import Counter
ROOT=Path(__file__).resolve().parent
RING=ROOT.parent/'sync/unfair-01-polynomials/computation/results/051-s12/frozen/old_ring.py'
s=importlib.util.spec_from_file_location('ring',RING);m=importlib.util.module_from_spec(s);s.loader.exec_module(m)
POINTS=((1,1),(1,2),(2,1),(1,4),(4,1),(3,4),(4,3),(6,5),(5,6))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def atomic(p,b):
    tmp=p.with_suffix(p.suffix+'.partial')
    with tmp.open('wb') as f:f.write(b);f.flush();os.fsync(f.fileno())
    tmp.replace(p)
def schur(a):
    a=list(a)
    while len(a)>1:
        div=math.gcd(*a)
        if a[-1]<0:div=-div
        a=[v//div for v in a]
        if abs(a[0])>=a[-1]:return len(a)-1
        n=len(a)-1;first,last=a[0],a[-1];a=[last*a[i]-first*a[n-i] for i in range(1,n+1)]
    return None
def norm(P,D):
    H=[(0,0)]*15;Q=[m.conj(z,D) for z in P]
    for i,x in enumerate(P):
        for j,y in enumerate(Q):H[i+j]=m.add(H[i+j],m.mul(x,y,D))
    assert all(z[1]==0 for z in H);return [z[0] for z in H]
def worker(inp,out,D,code,deadline):
    inp,out=Path(inp),Path(out);target=out/f'field-D{D}.json';source=inp/f'field-D{D}.json'
    while not source.exists() and time.time()<deadline:time.sleep(2)
    if not source.exists():return {'D':D,'status':'partial_missing_input'}
    binding={'source_sha256':code,'input_sha256':sha(source),'ring_sha256':sha(RING)}
    if target.exists():
        x=json.loads(target.read_text());assert x['binding']==binding and x['status']=='complete' and sha(out/x['ledger_file'])==x['ledger_sha256'];return {'D':D,'status':'complete_cached','counts':x['counts']}
    begin=time.monotonic();x=json.loads(source.read_text());assert x['status']=='complete' and x['D']==D
    mothers={}
    for row in x['rows']:
        key=tuple(row[:6]+[row[14]]);mask=mothers.get(key,0)
        for f in row[17]:assert -8<=f<=8;mask|=1<<(f+8)
        mothers[key]=mask
    del x;counts=Counter(unique_bases=len(mothers));codes=bytearray();sur=[]
    for (L,am,an,k,dm,dn,e),mask in sorted(mothers.items()):
        a=(am,an);d=(dm,dn);b=m.add(m.scale(L,a),(k,0))
        tests=[]
        for p,q in POINTS:
            w=[p**i*q**(7-i) for i in range(8)]
            z=m.add((w[0]+L*w[1]+e*w[5]+w[7],0),m.add(m.add(m.scale(w[2],a),m.scale(w[3],b)),m.scale(w[4],d)))
            tests.append((z,w[6]))
        for f in range(-8,9):
            if not mask&(1<<(f+8)):continue
            counts['unique_mothers']+=1;codebyte=None
            for j,(base,fc) in enumerate(tests):
                z=m.add(base,(fc*f,0))
                if m.sign(z,D)<=0 or m.sign(m.conj(z,D),D)<=0:codebyte=10+j;counts[f'axis{j}']+=1;break
            if codebyte is None:
                counts['positive_axis']+=1;H=norm([(1,0),(L,0),a,b,d,(e,0),(f,0),(1,0)],D)
                degree=schur([v<<i for i,v in enumerate(H)])
                if degree is not None:codebyte=20+degree;counts['outer_rejected']+=1
                else:
                    degree=schur([v<<i for i,v in enumerate(H[::-1])])
                    if degree is not None:codebyte=40+degree;counts['inner_rejected']+=1
                    else:codebyte=255;counts['survived']+=1;sur.append([D,L,am,an,k,dm,dn,e,f])
            codes.append(codebyte)
        if time.time()>deadline:
            # No partial byte sequence receives a complete marker.
            return {'D':D,'status':'partial_time_limit','processed':counts['unique_mothers']}
    assert len(codes)==counts['unique_mothers'];packed=gzip.compress(bytes(codes),mtime=0);ledger=f'field-D{D}-codes.gz';atomic(out/ledger,packed)
    result={'status':'complete','D':D,'binding':binding,'counts':dict(counts),'survivors':sur,'ledger_file':ledger,'ledger_sha256':hashlib.sha256(packed).hexdigest(),'raw_code_sha256':hashlib.sha256(codes).hexdigest(),'seconds':time.monotonic()-begin};atomic(target,(json.dumps(result,separators=(',',':'))+'\n').encode());return {'D':D,'status':'complete','counts':dict(counts),'seconds':result['seconds']}
def main():
    if not __debug__:raise RuntimeError('no -O')
    ap=argparse.ArgumentParser();ap.add_argument('--input',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);ap.add_argument('--workers',type=int,default=4);ap.add_argument('--fields',type=int,nargs='*');ap.add_argument('--seconds',type=int,default=14400);a=ap.parse_args();a.out.mkdir(exist_ok=True,parents=True)
    fields=json.loads((a.input/'summary.json').read_text())['fields'];fields=fields if a.fields is None else a.fields;code=sha(Path(__file__));deadline=time.time()+a.seconds;done=[]
    with cf.ProcessPoolExecutor(max_workers=a.workers) as pool:
        tasks=[pool.submit(worker,str(a.input),str(a.out),D,code,deadline) for D in fields]
        for future in cf.as_completed(tasks):
            z=future.result();done.append(z);print(json.dumps(z),flush=True)
            atomic(a.out/'progress.json',(json.dumps({'completed_fields':len(done),'expected_fields':len(fields),'records':done})+'\n').encode())
    counts=Counter()
    for z in done:counts.update(z.get('counts',{}))
    result={'status':'complete' if all(z['status'].startswith('complete') for z in done) else 'partial_unknown','fields':fields,'counts':dict(counts),'source_sha256':code,'scope':'Necessary positive-axis/root-annulus only; survivors require full cofactor trees.'};atomic(a.out/'summary.json',(json.dumps(result)+'\n').encode());print(json.dumps(result),flush=True)
if __name__=='__main__':main()
