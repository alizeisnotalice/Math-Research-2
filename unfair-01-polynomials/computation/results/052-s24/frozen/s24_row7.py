"""Complete c/e fifth/seventh-row necessary screen; f kept as a set.

No fixed shift window. Integer O_K arithmetic; per-field checkpoints.
Survivors still require full mother and arbitrary cofactor continuation.
"""
import argparse,hashlib,importlib.util,json,math,os,sys,time
from pathlib import Path
from collections import Counter,defaultdict
ROOT=Path(__file__).resolve().parent
RING=ROOT/'old_ring.py'
sp=importlib.util.spec_from_file_location('ring',RING);m=importlib.util.module_from_spec(sp);sp.loader.exec_module(m)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def floor_ring(x,D):
    R,V,den=m.rv(x,D);B=abs(R)+abs(V)*(D+1)+2;lo,hi=-B,B
    while hi-lo>1:
        mid=(lo+hi)//2
        if m.sign(m.add(x,(-mid,0)),D)>=0:lo=mid
        else:hi=mid
    return lo
def shifts(x,D):
    n=-floor_ring(x,D)
    return [j for j in (n-1,n,n+1) if unit(m.add(x,(j,0)),D)]
def unit(x,D):return m.sign(x,D)>=0 and m.sign(m.add((1,0),m.neg(x)),D)>=0
def binrow(U,V,cross,D):
    base=m.add(m.add(U,V),cross);assert base[1]==0
    out=[]
    for g in shifts(U,D):
        for C in (0,1):
            q=C-base[0]-g
            if unit(m.add(V,(q,0)),D):out.append((g,q))
    assert len(out)==len(set(out));return out
def one(rows):
    counts=Counter();out=[]
    for z in rows:
        D=z['D'];a=(z['a_m'],z['a_n']);d=(z['d_m'],z['d_n']);ab=m.conj(a,D);db=m.conj(d,D);g2,q2,g4,q4=[z[k] for k in ('g2','q2','g4','q4')];a2=m.mul(a,a,D);a3=m.mul(a2,a,D);ad=m.mul(a,d,D)
        A2=m.add(a,(g2,0));B2=m.add(ab,(q2,0))
        for c in range(-8,9):
            counts['raw_c']+=1
            fifth=binrow(m.scale(c,m.add((g2,0),m.neg(a))),m.scale(c,m.add((q2,0),m.neg(ab))),(0,0),D)
            if not fifth:counts['row5_rejected_c']+=1;continue
            counts['row5_states']+=len(fifth)
            # s6=-2a^3+6ad+3c^2-6f. e does not occur.
            s6base=m.add(m.add(m.scale(-2,a3),m.scale(6,ad)),(3*c*c,0))
            fs=[f for f in range(-8,9) if m.abs_both_lt(m.add(s6base,(-6*f,0)),D,12355,64)]
            if not fs:counts['no_f_c']+=1;continue
            for eta,zeta in fifth:
                A5=m.add(m.scale(c,m.add((g2,0),m.neg(a))),(eta,0));B5=m.add(m.scale(c,m.add((q2,0),m.neg(ab))),(zeta,0));cross=m.add(m.mul(A2,B5,D),m.mul(A5,B2,D))
                assert m.sign(cross,D)>=0
                for e in range(-29,30):
                    counts['raw_e_states']+=1
                    s5=m.add(m.scale(5*c,a),(-5*e,0))
                    if not m.abs_both_lt(s5,D,3107,32):counts['s5_rejected']+=1;continue
                    counts['s5_pass']+=1
                    s7=m.add(m.add(m.scale(-7*c,a2),m.scale(7*e,a)),m.add(m.scale(7*c,d),(-7,0)))
                    if not m.abs_both_lt(s7,D,49283,128):counts['s7_rejected']+=1;continue
                    counts['s7_pass']+=1
                    U=m.add((1+c*g4+e*g2,0),m.add(m.scale(eta-e,a),m.scale(-c,d)))
                    V=m.add((1+c*q4+e*q2,0),m.add(m.scale(zeta-e,ab),m.scale(-c,db)))
                    seventh=binrow(U,V,cross,D)
                    if not seventh:counts['row7_rejected']+=1;continue
                    counts['row7_states']+=len(seventh)
                    for g7,q7 in seventh:
                        out.append([z['a_m'],z['a_n'],z['d_m'],z['d_n'],g2,q2,g4,q4,z['rho6'],z['sigma6'],c,e,eta,zeta,g7,q7,fs]);counts['mothers_with_f']+=len(fs)
    return dict(counts),out
def atomic(p,x):
    temp=p.with_suffix('.partial')
    with temp.open('w') as f:json.dump(x,f,separators=(',',':'));f.write('\n');f.flush();os.fsync(f.fileno())
    temp.replace(p)
def main():
    if not __debug__:raise RuntimeError('no -O')
    ap=argparse.ArgumentParser();ap.add_argument('--input',type=Path,default=ROOT/'s24-row6.json');ap.add_argument('--out',type=Path,required=True);ap.add_argument('--seconds',type=float,default=14400);args=ap.parse_args();args.out.mkdir(parents=True,exist_ok=True);x=json.loads(args.input.read_text());groups=defaultdict(list)
    for z in x['rows']:groups[z['D']].append(z)
    assert sum(map(len,groups.values()))==20849
    binding={'input_sha256':sha(args.input),'source_sha256':sha(Path(__file__)),'ring_sha256':sha(RING)};start=time.monotonic();done=[];tot=Counter();number=0
    for D,rows in sorted(groups.items()):
        p=args.out/f'field-D{D}.json'
        if p.exists():v=json.loads(p.read_text());assert v['binding']==binding and v['D']==D
        else:
            if time.monotonic()-start>args.seconds:break
            counts,sur=one(rows);v={'status':'complete','D':D,'binding':binding,'counts':counts,'rows':sur};atomic(p,v)
        tot.update(v['counts']);number+=len(v['rows']);done.append(D);print(json.dumps({'D':D,'done_fields':len(done),'row7_states':number,'elapsed':time.monotonic()-start}),flush=True)
    o={'status':'complete' if len(done)==len(groups) else 'partial_unknown','binding':binding,'fields':done,'expected_fields':len(groups),'counts':dict(tot),'row7_states':number,'seconds':time.monotonic()-start,'python':sys.version,'scope':'Necessary c/e/f and rows<=7 only; no root annulus or full cofactor exclusion.'};atomic(args.out/'summary.json',o);print(json.dumps(o),flush=True)
if __name__=='__main__':main()
