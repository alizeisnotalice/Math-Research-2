#!/usr/bin/env python3
"""Independent exact necessary row-5/row-7 sieve for the S={2,4} model.

Input rows retain c in [-8,8] and f affinely.  This program enumerates c and
e in [-29,29], but never assigns f.  Thus every retained object is only a
necessary coefficient prefix, never a claimed factorisation or counterexample.
All integer translations are found with an unbounded exact quadratic floor.
"""
import argparse, hashlib, json, os, sys, time
from collections import Counter, defaultdict
from pathlib import Path

INPUT = Path('/Volumes/data/CodexTemp/T2/unfair-01-polynomials/agents/algebraic_global/round052-s24/row6-full/summary.json')
S5_BOUND = (3107, 32)
S6_BOUND = (12355, 64)
S7_BOUND = (49283, 128)

def add(x, y): return (x[0] + y[0], x[1] + y[1])
def neg(x): return (-x[0], -x[1])
def scl(k, x): return (k*x[0], k*x[1])
def mul(x, y, D):
    m,n=x; u,v=y
    if D % 4 == 1: return (m*u + n*v*((D-1)//4), m*v+n*u+n*v)
    return (m*u + D*n*v, m*v+n*u)
def bar(x, D): return (x[0]+x[1], -x[1]) if D % 4 == 1 else (x[0], -x[1])
def rv(x, D): return (2*x[0]+x[1], x[1], 2) if D % 4 == 1 else (x[0], x[1], 1)
def sign(x, D):
    R,V,s = rv(x,D)
    if V == 0: return (R > 0) - (R < 0)
    if V > 0:
        if R >= 0: return 1
        return 1 if V*V*D > R*R else -1
    if R <= 0: return -1
    return 1 if R*R > V*V*D else -1
def floorq(x, D):
    lo,hi=-1,1
    while sign(add(x,(-lo,0)),D) < 0: hi=lo; lo*=2
    while sign(add(x,(-hi,0)),D) >= 0: lo=hi; hi*=2
    while hi-lo > 1:
        mid=(lo+hi)//2
        if sign(add(x,(-mid,0)),D) >= 0: lo=mid
        else: hi=mid
    return lo
def shifts(base, upper, D):
    lo=-floorq(base,D)
    hi=floorq(add(upper,neg(base)),D)
    return list(range(lo,hi+1)) if lo <= hi else []
def abs_lt(x, frac, D):
    """Strict |x|<num/den at both real embeddings, using integers only."""
    num,den=frac
    R,V,s=rv(x,D)
    room=s*num-den*abs(R)
    return room > 0 and (den*V)*(den*V)*D < room*room
def atom(path, data):
    t=path.with_suffix(path.suffix+'.partial')
    with t.open('wb') as f: f.write(data); f.flush(); os.fsync(f.fileno())
    os.replace(t,path)

def c5_options(r, c):
    """Return (eta5,theta5,C5,A5,B5) exact choices for this c."""
    D=r['D']; a=(r['a_m'],r['a_n']); aa=bar(a,D)
    g2,q2=r['g2'],r['q2']
    if c == 0:
        return [(0,0,0,(0,0),(0,0)), (1,0,1,(1,0),(0,0)), (0,1,1,(0,0),(1,0))]
    ba=add((c*g2,0), neg(scl(c,a)))
    bb=add((c*q2,0), neg(scl(c,aa)))
    eta=shifts(ba,(1,0),D)
    if len(eta) > 1: raise AssertionError(('eta nonunique',r,c,ba,eta))
    if not eta: return []
    eta=eta[0]
    theta0=add(add((1,0),neg(ba)),neg(bb))
    theta0=add(theta0,(-eta,0))
    if theta0[1] != 0: raise AssertionError(('theta irrational',r,c,theta0))
    theta=theta0[0]
    A5=add(ba,(eta,0)); B5=add(bb,(theta,0))
    if not (sign(A5,D)>0 and sign(B5,D)>0 and add(A5,B5)==(1,0)):
        raise AssertionError(('strict row5',r,c,eta,theta,A5,B5))
    return [(eta,theta,1,A5,B5)]

def f_options(r, c):
    """All f in the finite root-bound box passing the two s6 bounds."""
    D=r['D']; a=(r['a_m'],r['a_n']); h=(r['h_m'],r['h_n'])
    fixed=add(add(scl(-2,mul(mul(a,a,D),a,D)),scl(6,mul(a,h,D))),(3*c*c,0))
    return [f for f in range(-8,9) if abs_lt(add(fixed,(-6*f,0)),S6_BOUND,D)]

def one(r, c, e, fvals):
    D=r['D']; a=(r['a_m'],r['a_n']); h=(r['h_m'],r['h_n'])
    aa,hh=bar(a,D),bar(h,D)
    g2,q2,g4,q4=r['g2'],r['q2'],r['g4'],r['q4']
    # Newton necessary bounds; both embeddings are checked by abs_lt.
    s5=add(scl(5*c,a),(-5*e,0))
    if not abs_lt(s5,S5_BOUND,D): return [], 's5'
    s7=add(add(add(scl(-7*c,mul(a,a,D)),scl(7*e,a)),scl(7*c,h)),(-7,0))
    if not abs_lt(s7,S7_BOUND,D): return [], 's7'
    out=[]
    for eta,theta,C5,A5,B5 in c5_options(r,c):
        cross=add(mul(add(a,(g2,0)),B5,D),mul(A5,add(aa,(q2,0)),D))
        baseA=add(add(add(add((1,0),scl(eta-e,a)),scl(c,(g4,0))),neg(scl(c,h))),scl(e,(g2,0)))
        baseB=add(add(add(add((1,0),scl(theta-e,aa)),scl(c,(q4,0))),neg(scl(c,hh))),scl(e,(q2,0)))
        for C7 in (0,1):
            target=add((C7,0),neg(cross))
            if sign(target,D) < 0: continue
            for g7 in shifts(baseA,target,D):
                q7v=add(add(target,neg(baseA)),neg(baseB)); q7v=add(q7v,(-g7,0))
                if q7v[1] != 0: raise AssertionError(('q7 irrational',r,c,e,eta,theta,C7,q7v))
                q7=q7v[0]; A7=add(baseA,(g7,0)); B7=add(baseB,(q7,0))
                if not (sign(A7,D) >= 0 and sign(B7,D) >= 0 and add(A7,B7)==target):
                    raise AssertionError(('row7',r,c,e,eta,theta,C7,g7,q7,A7,B7,target))
                z=dict(r)
                z.update({'c':c,'e':e,'eta5':eta,'theta5':theta,'C5':C5,
                          'G5_equals':'eta5-e','Q5_equals':'theta5-e',
                          'G6_equals':'R6-f+c^2','Q6_equals':'S6-f+c^2',
                          'f_values_passing_s6':fvals,
                          'g7':g7,'q7':q7,'C7':C7,
                          'scope_note':'necessary prefix through row 7; f remains free and this is not a factorization'})
                out.append(z)
    return out, 'prefix'

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--out',type=Path,required=True)
    ap.add_argument('--start',type=int,default=0); ap.add_argument('--count',type=int,default=None)
    ap.add_argument('--seconds',type=float,default=14400); z=ap.parse_args(); z.out.mkdir(parents=True,exist_ok=True)
    ibytes=INPUT.read_bytes(); inp=json.loads(ibytes)
    assert inp['status']=='complete' and len(inp['rows'])==20849
    groups=defaultdict(list)
    for r in inp['rows']: groups[r['D']].append(r)
    fs=sorted(groups); assert len(fs)==53
    if z.count is None: z.count=len(fs)
    assert 0 <= z.start <= len(fs) and 0 <= z.count and z.start+z.count <= len(fs)
    source=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    bind0={'source_sha256':source,'input_sha256':hashlib.sha256(ibytes).hexdigest()}
    begun=time.monotonic(); done=[]
    for i in range(z.start,z.start+z.count):
        D=fs[i]; path=z.out/f'field-{i:03d}-D{D}.json'; bind={**bind0,'index':i,'D':D,'bases':len(groups[D])}
        if path.exists():
            old=json.loads(path.read_text()); assert old['status']=='complete' and old['binding']==bind; done.append(i); continue
        if time.monotonic()-begun > z.seconds: break
        ct=Counter(); rows=[]; then=time.monotonic()
        for r in groups[D]:
            for c in range(-8,9):
                fvals=f_options(r,c)
                ct['raw_c'] += 1
                if not fvals:
                    ct['no_f_c'] += 1
                    continue
                for e in range(-29,30):
                    got,why=one(r,c,e,fvals); ct[why]+=1; rows.extend(got)
        obj={'status':'complete','binding':bind,'counts':dict(ct),'rows':rows,'elapsed_seconds':time.monotonic()-then}
        atom(path,(json.dumps(obj,sort_keys=True,separators=(',',':'))+'\n').encode()); done.append(i)
        print(json.dumps({'index':i,'D':D,'counts':dict(ct),'rows':len(rows),'seconds':obj['elapsed_seconds']},sort_keys=True),flush=True)
    total=Counter(); rows=[]
    for i in done:
        D=fs[i]; x=json.loads((z.out/f'field-{i:03d}-D{D}.json').read_text()); total.update(x['counts']); rows.extend(x['rows'])
    obj={'status':'complete' if len(done)==z.count else 'partial_unknown',
         'scope':'Exact necessary S={2,4} rows 5 and 7; c,e enumerated in root-bound boxes; f remains affine. Survivors are not polynomial factorizations.',
         'completed_indices':done,'counts':dict(total),'rows':rows,'binding':bind0,
         'elapsed_seconds':time.monotonic()-begun,'python':sys.version}
    atom(z.out/'summary.json',(json.dumps(obj,sort_keys=True,indent=2)+'\n').encode())
    print(json.dumps({'status':obj['status'],'fields':len(done),'counts':dict(total),'rows':len(rows),'seconds':obj['elapsed_seconds']},sort_keys=True))
if __name__ == '__main__': main()
