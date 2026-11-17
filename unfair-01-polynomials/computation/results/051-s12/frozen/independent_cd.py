#!/usr/bin/env python3
"""Independent exact c,d necessary screen from S12 exact-floor low rows.

No parent CD code is imported or read.  Input is only the frozen S12 low-row
ledger.  All comparisons use integer sqrt(D) sign tests; e and f are not
enumerated.  A5=B5=0 is represented at e=0 after proving its nonrational
part vanishes through the recorded high-end relation.
"""
from __future__ import annotations
import argparse, hashlib, json, math, os, sys, time
from collections import Counter, defaultdict
from pathlib import Path

INPUT=Path('/Volumes/data/CodexTemp/T2/unfair-01-polynomials/agents/algebraic_global/round051-s12/ab-lowprefix-floor-exact/full/summary.json')

def add(x,y):return(x[0]+y[0],x[1]+y[1])
def neg(x):return(-x[0],-x[1])
def scale(k,x):return(k*x[0],k*x[1])
def mul(x,y,D):return(x[0]*y[0]+x[1]*y[1]*D,x[0]*y[1]+x[1]*y[0])
def conjugate(x):return(x[0],-x[1])
def from_omega(m,n,D):return(2*m+n,n,2) if D%4==1 else (m,n,1)
def to_pair(m,n,D):
    R,V,s=from_omega(m,n,D)
    # Numerator pair representing (R+V sqrt(D))/s.  We retain denominator s
    # externally only here; all later values use denominator-1 scaled pairs.
    if s==1:return(R,V)
    # For D=1 mod4 we use doubled coordinates globally; adjust products below
    # by converting to rational-pair triples instead.  This function is unused.
    raise AssertionError

# General rational quadratic pair (r+s sqrt(D))/den, den>0, normalized only
# when convenient.  Fractions are avoided: every operation tracks denominator.
def qadd(x,y):
    r,s,d=x;u,v,e=y
    return(r*e+u*d,s*e+v*d,d*e)
def qneg(x):return(-x[0],-x[1],x[2])
def qscale(k,x):return(k*x[0],k*x[1],x[2])
def qmul(x,y,D):
    r,s,d=x;u,v,e=y
    return(r*u+s*v*D,r*v+s*u,d*e)
def qbar(x):return(x[0],-x[1],x[2])
def qint(k):return(k,0,1)
def qsign(x,D):
    r,s,d=x; assert d>0
    if s==0:return(r>0)-(r<0)
    if s>0:
        if r>=0:return 1
        return 1 if s*s*D>r*r else -1
    if r<=0:return-1
    return 1 if r*r>s*s*D else -1
def qeq(x,y):
    return x[0]*y[2]==y[0]*x[2] and x[1]*y[2]==y[1]*x[2]
def qisint(x):return x[1]==0 and x[0]%x[2]==0
def qinteger(x):assert qisint(x);return x[0]//x[2]
def qabsbothlt(x,D,num,den):
    # max(|r +/- s sqrtD|)/d < num/den
    r,s,d=x;room=d*num-den*abs(r)
    return room>0 and (den*s)*(den*s)*D<room*room
def qfloor(x,D):
    lo,hi=-1,1
    while qsign(qadd(x,qint(-lo)),D)<0:hi=lo;lo*=2
    while qsign(qadd(x,qint(-hi)),D)>=0:lo=hi;hi*=2
    while hi-lo>1:
        mid=(lo+hi)//2
        if qsign(qadd(x,qint(-mid)),D)>=0:lo=mid
        else:hi=mid
    return lo
def all_shifts(bases,D):
    """All integer z putting every base in [0,1]; exact and unbounded."""
    lo=max(-qfloor(x,D) for x in bases)
    hi=min(qfloor(qadd(qint(1),qneg(x)),D) for x in bases)
    if lo>hi:return []
    if hi-lo>1:raise AssertionError(('more than two endpoint shifts',bases,lo,hi))
    ans=list(range(lo,hi+1))
    for z in ans:
        for x in bases:
            v=qadd(x,qint(z));assert qsign(v,D)>=0 and qsign(qadd(qint(1),qneg(v)),D)>=0
    return ans
def omega(m,n,D):return(2*m+n,n,2) if D%4==1 else(m,n,1)

def one(row):
    D=row['D'];a=omega(row['a_m'],row['a_n'],D);b=omega(row['b_m'],row['b_n'],D)
    aa,bb=qbar(a),qbar(b)
    g1,q1,g2,q2=row['g1'],row['q1'],row['g2'],row['q2'];alpha,beta=row['alpha'],row['beta']
    c0=row['c_residue_mod_alpha'];out=[];unit_cases=[];counts=Counter(lowrow=1)
    for c in range(-52,53):
        counts['c_raw']+=1
        if c%alpha!=c0:continue
        counts['c_congruence']+=1
        g3=row['rho']-c;q3=row['q0']-c
        assert g3%alpha==0 and q3%alpha==0
        g4=-beta*(g3//alpha);q4=-beta*(q3//alpha)
        # s3=-a^3+3ab-3c, prior to choosing d.
        s3=qadd(qadd(qscale(-1,qmul(qmul(a,a,D),a,D)),qscale(3,qmul(a,b,D))),qint(-3*c))
        if not qabsbothlt(s3,D,203,8):counts['s3_rejected']+=1;continue
        counts['s3_pass']+=1
        # A4=d+baseA, B4=d+baseB.  The exact interval intersection preserves
        # both 0 and 1 endpoint shifts, hence may have two integers.
        baseA=qadd(qadd(qadd(qscale(g3,a),qscale(g2,b)),qint(c*g1)),qint(g4))
        baseB=qadd(qadd(qadd(qscale(q3,aa),qscale(q2,bb)),qint(c*q1)),qint(q4))
        ds=[z for z in all_shifts((baseA,baseB),D) if -52<=z<=52]
        counts['d_interval_candidates']+=len(ds)
        if not ds:counts['d_interval_empty']+=1;continue
        # Lower rows, needed for C4/C5.  These are recomputed, not trusted.
        A1=qadd(a,qint(g1));B1=qadd(aa,qint(q1))
        A2=qadd(qadd(b,qscale(g1,a)),qint(g2));B2=qadd(qadd(bb,qscale(q1,aa)),qint(q2))
        A3=qadd(qadd(qadd(qint(c),qscale(g1,b)),qscale(g2,a)),qint(g3))
        B3=qadd(qadd(qadd(qint(c),qscale(q1,bb)),qscale(q2,aa)),qint(q3))
        for d in ds:
            A4=qadd(baseA,qint(d));B4=qadd(baseB,qint(d))
            assert qsign(A4,D)>=0 and qsign(B4,D)>=0 and qsign(qadd(qint(1),qneg(A4)),D)>=0 and qsign(qadd(qint(1),qneg(B4)),D)>=0
            # s4=a^4-4a^2b+2b^2+4ac-4d.
            s4=qadd(qadd(qadd(qadd(qmul(qmul(qmul(a,a,D),a,D),a,D),qscale(-4,qmul(qmul(a,a,D),b,D))),qscale(2,qmul(b,b,D))),qscale(4*c,a)),qint(-4*d))
            C4=qadd(qadd(A4,B4),qadd(qadd(qmul(A1,B3,D),qmul(A2,B2,D)),qmul(A3,B1,D)))
            if not qisint(C4):raise AssertionError(('C4 not integral',row,c,d,C4))
            unit_cases.append({'D':D,'a_m':row['a_m'],'a_n':row['a_n'],'b_m':row['b_m'],'b_n':row['b_n'],'T':row['T'],'N':row['N'],'c':c,'d':d,
                               'g0':1,'g1':g1,'g2':g2,'g3':g3,'g4':g4,'q0':1,'q1':q1,'q2':q2,'q3':q3,'q4':q4,
                               'A4_pair':[A4[0],A4[1],A4[2]],'B4_pair':[B4[0],B4[1],B4[2]],'s4_pass':qabsbothlt(s4,D,787,16),'C4':qinteger(C4)})
            s4ok=qabsbothlt(s4,D,787,16)
            if s4ok:counts['s4_pass']+=1
            else:counts['s4_rejected']+=1
            if not (qisint(C4) and qinteger(C4) in (0,1)):
                # Kept independent of the Newton check, so this is a
                # readable all-unit witness count rather than a sequential
                # post-s4 count.
                counts['C4_rejected_all_unit']+=1;continue
            counts['C4_pass']+=1
            if not s4ok:counts['C4_pass_but_s4_rejected']+=1;continue
            C5=qadd(qadd(qadd(qmul(A1,B4,D),qmul(A2,B3,D)),qmul(A3,B2,D)),qmul(A4,B1,D))
            if not (qisint(C5) and qinteger(C5) in (0,1)):
                counts['C5_rejected']+=1;continue
            # e=0 representative: A5=B5=0 determines the signed integer rows.
            RA=qadd(qadd(qadd(qscale(g4,a),qscale(g3,b)),qint(c*g2+d*g1)),qint(0))
            RB=qadd(qadd(qadd(qscale(q4,aa),qscale(q3,bb)),qint(c*q2+d*q1)),qint(0))
            if not(qisint(RA) and qisint(RB)):raise AssertionError(('A5 projection failure',row,c,d,RA,RB))
            g5,q5=-qinteger(RA),-qinteger(RB)
            out.append({'D':D,'a_m':row['a_m'],'a_n':row['a_n'],'b_m':row['b_m'],'b_n':row['b_n'],'T':row['T'],'N':row['N'],'c':c,'d':d,
                        'g0':1,'g1':g1,'g2':g2,'g3':g3,'g4':g4,'g5_e0':g5,
                        'q0':1,'q1':q1,'q2':q2,'q3':q3,'q4':q4,'q5_e0':q5,
                        'C4':qinteger(C4),'C5':qinteger(C5)})
            counts['survivor']+=1
    return dict(counts),out,unit_cases

def atomic(path,data):
    tmp=path.with_suffix(path.suffix+'.partial')
    with tmp.open('wb') as f:f.write(data);f.flush();os.fsync(f.fileno())
    os.replace(tmp,path)

def main():
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);p.add_argument('--start',type=int,default=0);p.add_argument('--count',type=int,default=None);p.add_argument('--seconds',type=float,default=14400);args=p.parse_args();args.out.mkdir(parents=True,exist_ok=True)
    ibytes=INPUT.read_bytes();inp=json.loads(ibytes);assert inp['status']=='complete' and len(inp['rows'])==7245
    groups=defaultdict(list)
    for r in inp['rows']:groups[r['D']].append(r)
    # Some of the 81 quadratic fields have no retained low row.  The input
    # summary itself is the complete 7,245-row universe for this stage.
    fields=sorted(groups);assert sum(len(groups[D]) for D in fields)==7245
    if args.count is None:args.count=len(fields)
    assert 0<=args.start and 0<=args.count and args.start+args.count<=len(fields)
    source=hashlib.sha256(Path(__file__).read_bytes()).hexdigest();binding={'source_sha256':source,'input_sha256':hashlib.sha256(ibytes).hexdigest()};start=time.monotonic();done=[]
    for i in range(args.start,args.start+args.count):
        D=fields[i];path=args.out/f'field-{i:02d}-D{D}.json';b={**binding,'index':i,'D':D,'low_rows':len(groups[D])}
        if path.exists():
            x=json.loads(path.read_text());assert x['status']=='complete' and x['binding']==b;done.append(i);continue
        if time.monotonic()-start>args.seconds:break
        ct=Counter();rows=[];unit_cases=[];begun=time.monotonic()
        for r in groups[D]:c,s,u=one(r);ct.update(c);rows.extend(s);unit_cases.extend(u)
        obj={'status':'complete','binding':b,'counts':dict(ct),'rows':rows,'unit_cases':unit_cases,'elapsed_seconds':time.monotonic()-begun};atomic(path,(json.dumps(obj,sort_keys=True,separators=(',',':'))+'\n').encode());done.append(i)
        print(json.dumps({'index':i,'D':D,'counts':dict(ct),'rows':len(rows),'unit_cases':len(unit_cases),'seconds':obj['elapsed_seconds']},sort_keys=True),flush=True)
    total=Counter();rows=[];unit_cases=[]
    for i in done:
        D=fields[i];x=json.loads((args.out/f'field-{i:02d}-D{D}.json').read_text());total.update(x['counts']);rows.extend(x['rows']);unit_cases.extend(x['unit_cases'])
    result={'status':'complete' if len(done)==args.count else 'partial_unknown','scope':'Necessary S12 c,d rows only; e=0 denotes an A5/B5 representative and no e/f completion is claimed.','completed_indices':done,'counts':dict(total),'rows':rows,'unit_cases':unit_cases,'binding':binding,'elapsed_seconds':time.monotonic()-start,'python':sys.version}
    atomic(args.out/'summary.json',(json.dumps(result,sort_keys=True,indent=2)+'\n').encode());print(json.dumps({'status':result['status'],'fields':len(done),'counts':dict(total),'rows':len(rows),'seconds':result['elapsed_seconds']},sort_keys=True))
if __name__=='__main__':main()
