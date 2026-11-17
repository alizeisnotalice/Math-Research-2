#!/usr/bin/env python3
"""Exact S={1,2} lattice/low-prefix screen, without the c,d,e,f mother box.

All elements are m+n*omega in the full ring of integers O_Q(sqrt(D)).  This
is deliberately not a Z[a] representation.  The screen records only
necessary conditions through C3 and the high-end A5=B5 projection relation.
It is not an existence search for G,Q and it never declares global exclusion.
"""
from __future__ import annotations

import argparse, hashlib, json, math, os, sys, time
from collections import Counter
from pathlib import Path


def squarefree(d: int) -> bool:
    p = 2
    while p * p <= d:
        if d % (p * p) == 0:
            return False
        p += 1
    return True


def mul(x, y, d):
    """Multiply m+n omega coordinates in the full integral basis."""
    m, n = x; u, v = y
    if d % 4 == 1:
        # omega^2=omega+(d-1)/4.
        return (m*u + n*v*((d-1)//4), m*v+n*u+n*v)
    return (m*u+d*n*v, m*v+n*u)


def add(x, y): return (x[0]+y[0], x[1]+y[1])
def neg(x): return (-x[0],-x[1])
def scale(k, x): return (k*x[0],k*x[1])


def conj(x, d):
    m, n = x
    return (m+n, -n) if d % 4 == 1 else (m, -n)


def trace(x, d):
    return 2*x[0]+x[1] if d % 4 == 1 else 2*x[0]


def norm(x, d):
    m, n = x
    return m*m+m*n+n*n*(1-d)//4 if d % 4 == 1 else m*m-d*n*n


def rv(x, d):
    m, n=x
    return (2*m+n,n,2) if d%4==1 else (m,n,1)


def sign(x, d):
    """Sign under the chosen embedding sqrt(D)>0; zero allowed."""
    R,V,s=rv(x,d)
    if V == 0:
        return (R>0)-(R<0)
    if V > 0:
        if R >= 0: return 1
        return 1 if V*V*d > R*R else -1
    if R <= 0: return -1
    return 1 if R*R > V*V*d else -1


def abs_both_lt(x, d, num, den):
    """Both embeddings have absolute value strictly below num/den."""
    R,V,s=rv(x,d)
    room=s*num-den*abs(R)
    return room>0 and (den*V)*(den*V)*d < room*room


def ring_candidates(d, num, den, positive_n):
    # The loop bounds are loose consequences of max embedding absolute value.
    if d%4==1:
        nb=(2*num+den-1)//den; mb=nb
    else:
        nb=(num+den-1)//den; mb=nb
    out=[]
    for n in range(1 if positive_n else -nb,nb+1):
        if n==0: continue
        for m in range(-mb,mb+1):
            x=(m,n)
            if abs_both_lt(x,d,num,den): out.append(x)
    return out


def field_list():
    out=[]
    for d in range(2,289):
        if not squarefree(d): continue
        aa=ring_candidates(d,17,2,True)
        if aa: out.append((d,aa,ring_candidates(d,117,4,False)))
    assert len(out)==81 and sum(len(a)*len(b) for _,a,b in out)==396738
    return out


def integer_shifts(base, upper, d, low=-110, high=110):
    """All z in the supplied safe finite range with 0<=base+z<=upper."""
    ans=[]
    for z in range(low,high+1):
        val=add(base,(z,0))
        if sign(val,d)>=0 and sign(add(upper,neg(val)),d)>=0: ans.append(z)
    if len(ans)>1: raise AssertionError(('interval too long',base,upper,ans))
    return ans


def check_pair(d,a,b):
    """Return first-stage code and, on full low/high survival, one ledger row."""
    # Newton s2=a^2-2b, evaluated in both embeddings.
    if not abs_both_lt(add(mul(a,a,d),scale(-2,b)),d,55,4): return 's2',None
    T,N=trace(a,d),norm(a,d)
    g1=-math.floor(0)  # overwritten below; makes its calculation visibly separate.
    # a lies in (-17/2,17/2).  Exact signs avoid converting sqrt(D) to float.
    for k in range(-9,10):
        if sign(add(a,(-k,0)),d)>=0 and sign(add(a,(-(k+1),0)),d)<0:
            g1=-k; break
    else: raise AssertionError(('floor failure',d,a))
    q1=1-T-g1
    A1=add(a,(g1,0)); B1=add(conj(a,d),(q1,0))
    if not (sign(A1,d)>0 and sign(add((1,0),neg(A1)),d)>0 and B1==add((1,0),neg(A1))):
        raise AssertionError(('row1',d,a,g1,q1,A1,B1))
    # S=1-A1*B1.  It is in (3/4,1), strictly since A1 is nonrational.
    S=add((1,0),neg(mul(A1,B1,d)))
    if not (sign(add(S,(-1,0)),d)<0 and sign(add(scale(4,S),(-3,0)),d)>0):
        raise AssertionError(('S bounds',d,a,S))
    A2base=add(b,scale(g1,a))
    shifts=integer_shifts(A2base,S,d)
    if not shifts: return 'A2',None
    g2=shifts[0]
    # Exact rational C2 identity: c2=Tr(b)+T*q1+g1+N-g1^2+g2+q2=1.
    q2=1-trace(b,d)-T*q1-g1-N+g1*g1-g2
    A2=add(A2base,(g2,0))
    B2=add(add(conj(b,d),scale(q1,conj(a,d))), (q2,0))
    if B2 != add(S,neg(A2)) or sign(A2,d)<0 or sign(B2,d)<0:
        raise AssertionError(('row2',d,a,b,g1,q1,g2,q2,A2,B2,S))
    K=add(mul(A1,B2,d),mul(A2,B1,d))
    if not (sign(K,d)>0 and sign(add((1,0),neg(K)),d)>0):
        raise AssertionError(('cross',d,K))
    # C3 without g3+q3 is 2c+Y.  The c coefficient is rational, so its
    # nonrational projection vanishes iff that of Y vanishes.
    Y=add(add(add(scale(g1,b),scale(q1,conj(b,d))),scale(g2,a)),scale(q2,conj(a,d)))
    Y=add(Y,K)
    if Y[1] != 0: return 'C3_projection',None
    Yint=Y[0]
    rho_shifts=integer_shifts(add(scale(g1,b),scale(g2,a)),add((1,0),neg(K)),d)
    if not rho_shifts: return 'A3',None
    rho=rho_shifts[0]                 # rho=c+g3, the unique low-row shift.
    q0=1-Yint-rho                    # q3=q0-c when C3=1.
    h=math.gcd(a[1],abs(b[1])); alpha=a[1]//h; beta=b[1]//h
    if (rho-q0)%alpha: return 'high_congruence',None
    row={
        'D':d,'a_m':a[0],'a_n':a[1],'b_m':b[0],'b_n':b[1],
        'T':T,'N':N,'g1':g1,'q1':q1,'g2':g2,'q2':q2,
        'alpha':alpha,'beta':beta,'rho':rho,'q0':q0,
        'c_residue_mod_alpha':rho%alpha,
    }
    return 'survivor',row


def atomic(path,data):
    temp=path.with_suffix(path.suffix+'.partial')
    with temp.open('wb') as f:
        f.write(data);f.flush();os.fsync(f.fileno())
    os.replace(temp,path)


def one_field(d,aa,bb):
    counts=Counter(raw=0); rows=[]
    for a in aa:
        for b in bb:
            counts['raw'] += 1
            code,row=check_pair(d,a,b); counts[code]+=1
            if row is not None: rows.append(row)
    assert counts['raw']==len(aa)*len(bb)
    assert sum(v for k,v in counts.items() if k!='raw')==counts['raw']
    return dict(counts),rows


def main():
    if not __debug__: raise RuntimeError('do not run with -O')
    ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True)
    ap.add_argument('--start',type=int,default=0);ap.add_argument('--count',type=int,default=81)
    ap.add_argument('--seconds',type=float,default=14400);args=ap.parse_args();args.out.mkdir(parents=True,exist_ok=True)
    fields=field_list();assert 0<=args.start<=len(fields) and 0<=args.count<=len(fields)-args.start
    source=hashlib.sha256(Path(__file__).read_bytes()).hexdigest();started=time.monotonic();done=[]
    for i in range(args.start,args.start+args.count):
        d,aa,bb=fields[i];path=args.out/f'field-{i:02d}-D{d}.json'
        binding={'source_sha256':source,'field_index':i,'D':d,'a_count':len(aa),'b_count':len(bb)}
        if path.exists():
            old=json.loads(path.read_text())
            if old.get('status')=='complete' and old.get('binding')==binding: done.append(i);continue
            raise RuntimeError(f'mismatched checkpoint {path}')
        if time.monotonic()-started>args.seconds: break
        begun=time.monotonic();counts,rows=one_field(d,aa,bb)
        out={'status':'complete','binding':binding,'counts':counts,'rows':rows,'elapsed_seconds':time.monotonic()-begun}
        atomic(path,(json.dumps(out,sort_keys=True,separators=(',',':'))+'\n').encode());done.append(i)
        print(json.dumps({'index':i,'D':d,'counts':counts,'rows':len(rows),'seconds':out['elapsed_seconds']},sort_keys=True),flush=True)
    total=Counter();rows=[]
    for i in done:
        d=fields[i][0];x=json.loads((args.out/f'field-{i:02d}-D{d}.json').read_text());total.update(x['counts']);rows.extend(x['rows'])
    result={'status':'complete' if len(done)==args.count else 'partial_unknown','scope':'S={1,2} O_K lattice pairs through C3 and A5=B5 projection only; no c,d,e,f box and no cofactor completion.',
            'requested_start':args.start,'requested_count':args.count,'completed_indices':done,'counts':dict(total),'rows':rows,
            'source_sha256':source,'python':sys.version,'elapsed_seconds':time.monotonic()-started,
            'codes':'s2: Newton s2 fails; A2: no second-row integer shift; C3_projection: Y not rational; A3: no third-row shift; high_congruence: no c residue satisfies both A5/B5 projection equations; survivor: all recorded necessary conditions.'}
    atomic(args.out/'summary.json',(json.dumps(result,sort_keys=True,indent=2)+'\n').encode())
    print(json.dumps({'status':result['status'],'completed':len(done),'counts':dict(total),'rows':len(rows),'seconds':result['elapsed_seconds']},sort_keys=True))


if __name__=='__main__': main()
