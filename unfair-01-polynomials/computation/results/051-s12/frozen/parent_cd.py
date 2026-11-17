"""Exact C4/C5 and Newton reduction after complete S12 low prefixes.

Tail e/f are NOT enumerated. e=0 here is only a representative for the
fifth integer cofactor coefficient, not removal from the full mother.
"""
import importlib.util,json,hashlib,time
from pathlib import Path
from collections import Counter
T=Path('/Volumes/data/CodexTemp/T2/unfair-01-polynomials')
SOURCE=T/'agents/algebraic_global/round051-s12/s12_ab_lowprefix_screen.py'
INPUT=T/'agents/algebraic_global/round051-s12/ab-lowprefix-floor-exact/full/summary.json'
sp=importlib.util.spec_from_file_location('ring',SOURCE);m=importlib.util.module_from_spec(sp);sp.loader.exec_module(m)
def floor_ring(x,D):
    R,V,den=m.rv(x,D);B=abs(R)+abs(V)*(D+1)+2;lo,hi=-B,B
    while hi-lo>1:
        mid=(lo+hi)//2
        if m.sign(m.add(x,(-mid,0)),D)>=0:lo=mid
        else:hi=mid
    return lo
def shifts(x,D):
    # All integers z making x+z in [0,1]. Rational endpoints both retained.
    n=-floor_ring(x,D)
    return [j for j in (n-1,n,n+1) if m.sign(m.add(x,(j,0)),D)>=0 and m.sign(m.add((1,0),m.neg(m.add(x,(j,0)))),D)>=0]
def conv(A,B,D,k):
    z=(0,0)
    for j in range(k+1):
        if j<len(A) and k-j<len(B):z=m.add(z,m.mul(A[j],B[k-j],D))
    return z
def powers(P,D,k):
    S=[(0,0)]
    for n in range(1,k+1):
        v=m.scale(n,P[n])
        for j in range(1,n):v=m.add(v,m.mul(P[j],S[n-j],D))
        S.append(m.neg(v))
    return S
def main():
    begin=time.monotonic();source=SOURCE.read_bytes();raw=INPUT.read_bytes();data=json.loads(raw);counts=Counter();rows=[];last194=[]
    for z in data['rows']:
        D=z['D'];a=(z['a_m'],z['a_n']);b=(z['b_m'],z['b_n']);alpha,beta=z['alpha'],z['beta']
        for c in range(-52,53):
            counts['c_raw']+=1
            if (c-z['rho'])%alpha or (c-z['q0'])%alpha:continue
            counts['c_congruence']+=1;g3=z['rho']-c;q3=z['q0']-c;g4=-beta*g3//alpha;q4=-beta*q3//alpha
            P=[(1,0),a,b,(c,0),(0,0)];S=powers(P,D,3)
            if not m.abs_both_lt(S[3],D,203,8):continue
            counts['s3']+=1;G=[1,z['g1'],z['g2'],g3,g4];Q=[1,z['q1'],z['q2'],q3,q4]
            base=m.add((g4+c*G[1],0),m.add(m.scale(g3,a),m.scale(G[2],b)))
            for d in shifts(base,D):
                counts['d_unit_candidates']+=1
                if not -52<=d<=52:continue
                counts['d_in_box']+=1
                P=[(1,0),a,b,(c,0),(d,0)];bar=[m.conj(x,D) for x in P]
                A=[conv(P,[(v,0) for v in G],D,k) for k in range(5)];B=[conv(bar,[(v,0) for v in Q],D,k) for k in range(5)]
                if any(m.sign(x,D)<0 or m.sign(m.add((1,0),m.neg(x)),D)<0 for x in A+B):counts['unit_rejected']+=1;continue
                counts['unit_passed']+=1
                last194.append({'prefix':z,'c':c,'d':d,'G':G,'Q':Q,'A':A,'B':B,'C4':conv(A,B,D,4)})
                if conv(A,B,D,4) not in ((0,0),(1,0)):counts['C4_rejected']+=1;continue
                counts['C4']+=1
                s4=powers(P,D,4)[4]
                if not m.abs_both_lt(s4,D,787,16):continue
                counts['s4']+=1
                # A5=B5=0; the omega part cancels by the high-end relation.
                x=m.add(m.scale(g4,a),m.add(m.scale(g3,b),(c*G[2]+d*G[1],0)))
                y=m.add(m.scale(q4,m.conj(a,D)),m.add(m.scale(q3,m.conj(b,D)),(c*Q[2]+d*Q[1],0)))
                assert x[1]==y[1]==0;G5=-x[0];Q5=-y[0]
                if conv(A+[(0,0)],B+[(0,0)],D,5) not in ((0,0),(1,0)):continue
                counts['C5']+=1;rows.append({'prefix':z,'c':c,'d':d,'G_at_e0':G+[G5],'Q_at_e0':Q+[Q5]})
    out={'status':'complete_necessary_cd_screen','scope':'Necessary prefix test only; no standalone claim of arbitrary cofactor exclusion until low-prefix coverage and formulas independently verified.','counts':dict(counts),'rows':rows,'last_unit_prefixes':last194,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'ring_source_sha256':hashlib.sha256(source).hexdigest(),'input_sha256':hashlib.sha256(raw).hexdigest(),'seconds':time.monotonic()-begin}
    Path(__file__).with_name('s12-cd-result.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items() if k not in ('rows','last_unit_prefixes')}))
if __name__=='__main__':main()
