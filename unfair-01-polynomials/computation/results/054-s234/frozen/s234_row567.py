"""Necessary rows5/6/7 and full e/f Newton box, exact integer O_K arithmetic.

Free G5/G6/G7 integer terms are represented after removing known integer
translations. No fixed window, no cofactor degree bound, no success theorem.
"""
import argparse,concurrent.futures as cf,functools,hashlib,importlib.util,json,os,sys,time
from pathlib import Path
from collections import Counter
ROOT=Path(__file__).resolve().parent
s=importlib.util.spec_from_file_location('r4',ROOT/'s234_row4.py');r4=importlib.util.module_from_spec(s);s.loader.exec_module(r4);m=r4.m
def binrow(U,V,cross,D,targets=(0,1)):
    base=m.add(m.add(U,V),cross);assert base[1]==0;out=[]
    for g in r4.shifts(U,D):
        for C in targets:
            q=C-base[0]-g
            if r4.unit(m.add(V,(q,0)),D):out.append((g,q))
    assert len(out)==len(set(out));return out
@functools.lru_cache(maxsize=65536)
def even(D,am,an,dm,dn,u2,v2,u4,v4):
    a=(am,an);d=(dm,dn);ab=m.conj(a,D);db=m.conj(d,D);A2=m.add(a,(u2,0));B2=m.add(ab,(v2,0));A4=m.add(m.add(d,m.scale(u2,a)),(u4,0));B4=m.add(m.add(db,m.scale(v2,ab)),(v4,0));cross=m.add(m.mul(A2,B4,D),m.mul(A4,B2,D))
    assert m.sign(cross,D)>0 and m.sign(m.add((1,0),m.neg(cross)),D)>0
    a2=m.mul(a,a,D);return a,ab,d,db,A2,B2,cross,a2,m.mul(a2,a,D),m.mul(a,d,D)
def one(D,source,target,code):
    source,target=Path(source),Path(target);inputhash=hashlib.sha256(source.read_bytes()).hexdigest();binding={'source_sha256':code,'input_sha256':inputhash,'ring_sha256':hashlib.sha256(r4.RING.read_bytes()).hexdigest(),'row4_code_sha256':hashlib.sha256((ROOT/'s234_row4.py').read_bytes()).hexdigest()}
    if target.exists():x=json.loads(target.read_text());assert x['binding']==binding and x['status']=='complete';return {'D':D,'counts':x['counts'],'cached':True}
    x=json.loads(source.read_text());assert x['status']=='complete' and x['D']==D;counts=Counter();out=[];start=time.monotonic()
    for L,am,an,k,dm,dn,u2,v2,u4,v4 in x['rows']:
        counts['row4_input']+=1;a=(am,an);d=(dm,dn)
        base5=m.add(m.add(m.scale(5*k,a),m.scale(5*L,d)),(-L**5-5*L*L*k,0))
        if not m.abs_both_lt(base5,D,7747,32):counts['s5_triangle_rejected']+=1;continue
        es=[e for e in range(-29,30) if m.abs_both_lt(m.add(base5,(-5*e,0)),D,3107,32)]
        if not es:counts['no_e_s5']+=1;continue
        counts['e_s5']+=len(es)
        a,ab,d,db,A2,B2,cross6,a2,a3,ad=even(D,am,an,dm,dn,u2,v2,u4,v4)
        U6=m.add(m.scale(u4+L*k,a),m.scale(u2+L*L,d));V6=m.add(m.scale(v4+L*k,ab),m.scale(v2+L*L,db));six=binrow(U6,V6,cross6,D,(1,));assert len(six)<=1
        if not six:counts['row6_rejected']+=1;continue
        rho,sigma=six[0];counts['row6_pass']+=1
        U5=m.add(m.scale(-k,a),m.scale(-L,d));V5=m.conj(U5,D);fifth=binrow(U5,V5,(0,0),D);assert len(fifth) in (1,3)
        rawg2=L*L+u2;rawq2=L*L+v2
        s6base=m.add(m.add(m.scale(-2,a3),m.scale(6,ad)),m.add(m.scale(-6*L*k,a),m.add(m.scale(-6*L*L,d),(L**6+6*L**3*k+3*k*k,0))))
        s7base=m.add(m.add(m.scale(-7*k,a2),m.scale(7*L*L*k,a)),m.add(m.scale(-7*L,ad),m.add(m.scale(7*(L**3+k),d),(-L**7-7*L**4*k-7*L*k*k-7,0))))
        for e in es:
            counts['e_before_f']+=1;s6e=m.add(s6base,(6*L*e,0));s7e=m.add(s7base,m.add(m.scale(7*e,a),(-7*L*L*e,0)))
            if not m.abs_both_lt(s6e,D,15427,64) or not m.abs_both_lt(s7e,D,49283+128*56*abs(L),128):counts['f_triangle_rejected']+=1;continue
            fs=[f for f in range(-8,9) if m.abs_both_lt(m.add(s6e,(-6*f,0)),D,12355,64) and m.abs_both_lt(m.add(s7e,(7*L*f,0)),D,49283,128)]
            if not fs:counts['no_f_s67']+=1;continue
            counts['e_with_f']+=1
            for eta,zeta in fifth:
                A5=m.add(U5,(eta,0));B5=m.add(V5,(zeta,0));cross7=m.add(m.mul(A2,B5,D),m.mul(A5,B2,D));U7=m.add(m.scale(eta-e-k*rawg2,a),m.scale(-L*rawg2-k,d));V7=m.add(m.scale(zeta-e-k*rawq2,ab),m.scale(-L*rawq2-k,db));seventh=binrow(U7,V7,cross7,D)
                if not seventh:counts['row7_rejected']+=1;continue
                for gamma,delta in seventh:
                    counts['row7_states']+=1;counts['state_f_combinations']+=len(fs);out.append([L,am,an,k,dm,dn,u2,v2,u4,v4,rho,sigma,eta,zeta,e,gamma,delta,fs])
    obj={'status':'complete','D':D,'binding':binding,'counts':dict(counts),'rows':out,'seconds':time.monotonic()-start};atomic(target,obj);even.cache_clear();return {'D':D,'counts':dict(counts),'seconds':obj['seconds']}
def atomic(p,x):
    tmp=p.with_suffix('.partial')
    with tmp.open('w') as f:json.dump(x,f,separators=(',',':'));f.write('\n');f.flush();os.fsync(f.fileno())
    tmp.replace(p)
def main():
    if not __debug__:raise RuntimeError('no -O')
    ap=argparse.ArgumentParser();ap.add_argument('--input',type=Path,default=ROOT/'row4-full');ap.add_argument('--out',type=Path,required=True);ap.add_argument('--workers',type=int,default=2);args=ap.parse_args();args.out.mkdir(parents=True,exist_ok=True);fields=json.loads((args.input/'summary.json').read_text())['fields'];assert len(fields)==58;code=hashlib.sha256(Path(__file__).read_bytes()).hexdigest();start=time.monotonic();done=[];counts=Counter()
    with cf.ProcessPoolExecutor(max_workers=args.workers) as pool:
        jobs=[pool.submit(one,D,str(args.input/f'field-D{D}.json'),str(args.out/f'field-D{D}.json'),code) for D in fields]
        for job in cf.as_completed(jobs):z=job.result();done.append(z['D']);counts.update(z['counts']);print(json.dumps(z),flush=True);atomic(args.out/'progress.json',{'done':done,'counts':dict(counts),'expected':58})
    assert sorted(done)==fields and counts['row4_input']==4459150
    result={'status':'complete','fields':fields,'counts':dict(counts),'source_sha256':code,'seconds':time.monotonic()-start,'python':sys.version,'record_schema':['ell','a_m','a_n','k','d_m','d_n','u2','v2','u4','v4','rho6','sigma6','eta5','zeta5','e','gamma7','delta7','f_values'],'scope':'Necessary complete e/f Newton and rows<=7 only; no root-annulus or arbitrary cofactor exclusion.'};atomic(args.out/'summary.json',result);print(json.dumps(result),flush=True)
if __name__=='__main__':main()
