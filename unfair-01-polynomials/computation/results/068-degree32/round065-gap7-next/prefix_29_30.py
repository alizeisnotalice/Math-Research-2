"""Guarded exact prefix enumeration for degrees 29 and 30."""
from collections import Counter
from dataclasses import dataclass
from itertools import product
from pathlib import Path
import json,sys

sys.dont_write_bytecode=True
from verify_degrees_25_26 import (T,ZERO,INTERVAL,setup,positive_low_rows,
                                  proof_exists_positive,linear_test,local_linear_rejection)
from verify_degree_19 import cross,element

CODES='ABXZ'
PARTNER={'A':'A','B':'B','X':'Z','Z':'F'}


def initial_shapes(j):
    for states in product('ABX','ABX','ABX',CODES,CODES,CODES,CODES):
        early=''.join(states)
        for mid in ('AB' if j%2==0 else (None,)):
            yield early,mid


def required_positive(early):
    a,b=positive_low_rows()
    for k,state in enumerate(early,8):
        if state in ('A','X'):a.add(k)
        if state in ('B','X'):b.add(k)
    return a,b


def masks(j,early):
    a,b=required_positive(early)
    return [((1,) if proof_exists_positive(a,b,k) else (0,1)) for k in range(12,j-7)]


def sequence(j,shape):
    early,mid,_=shape
    a=early+(mid or '')+''.join(PARTNER[z] for z in early[::-1])
    assert len(a)==j-15
    return a


def slots(j,shape):
    return tuple(k for k,z in enumerate(sequence(j,shape),8) if z in ('X','F'))


def model(t,j,shape,free):
    values=iter(free)
    residual=[]
    A,B=setup(t)
    for k,(side,c) in enumerate(zip(sequence(j,shape),[1]*4+list(shape[2]),strict=True),8):
        x=c-cross(A,B,k)
        if side=='A':a,b=x,0
        elif side=='B':a,b=0,x
        elif side in ('X','F'):
            a=next(values);b=x-a
        else:
            assert side=='Z'
            a=b=0
            residual.append((k,x))
        A.append(a);B.append(b)
    assert not list(values)
    return residual,A,B


@dataclass(frozen=True)
class Polynomial:
    coef:dict

    def __add__(self,r):
        r=polynomial(r)
        z=self.coef.copy()
        for a,b in r.coef.items():
            z[a]=z.get(a,ZERO)+b
            if z[a]==ZERO:
                z.pop(a)
        return Polynomial(z)

    __radd__=__add__

    def __neg__(self):
        return Polynomial({a:-b for a,b in self.coef.items()})

    def __sub__(self,r):return self+-polynomial(r)

    def __rsub__(self,r):return polynomial(r)+-self

    def __mul__(self,r):
        r=polynomial(r)
        z={}
        for a,b in self.coef.items():
            for c,d in r.coef.items():
                k=tuple(sorted(a+c))
                assert len(k)<=2
                z[k]=z.get(k,ZERO)+b*d
        return Polynomial({a:b for a,b in z.items() if b!=ZERO})

    __rmul__=__mul__


def polynomial(x):
    if isinstance(x,Polynomial):return x
    e=element(x)
    return Polynomial({():e} if e!=ZERO else {})


def exact_system(j,shape):
    n=len(slots(j,shape))
    free=[Polynomial({(i,):element(1)}) for i in range(n)]
    rows,_,_=model(polynomial(T),j,shape,free)
    assert len(rows)==n
    return [(k,polynomial(v).coef) for k,v in rows]


SCALE=1<<42


@dataclass(frozen=True)
class Box:
    lo:int
    hi:int

    def __add__(self,y):
        y=point(y)
        return Box(self.lo+y.lo,self.hi+y.hi)

    __radd__=__add__

    def __neg__(self):
        return Box(-self.hi,-self.lo)

    def __sub__(self,y):return self+-point(y)

    def __rsub__(self,y):return point(y)+-self

    def __mul__(self,y):
        y=point(y)
        v=(self.lo*y.lo,self.lo*y.hi,self.hi*y.lo,self.hi*y.hi)
        return Box(min(v)//SCALE,-(-max(v)//SCALE))

    __rmul__=__mul__


def point(x):
    if isinstance(x,Box):return x
    if isinstance(x,int):return Box(x*SCALE,x*SCALE)
    raise TypeError(type(x))


def canon(s):
    return s[0],s[1],tuple(s[2])


def t_box():
    lower=71942965*SCALE//10**8
    upper=-(-71942966*SCALE//10**8)
    return Box(lower,upper)


B0=point(0)


def fill_interval(S,state,strict=False):
    if state=='Z':
        return (B0,B0) if S.lo<=0<=S.hi else None
    limit=2*SCALE if state in ('X','F') else SCALE
    if S.hi<0 or S.lo>limit or strict and S.hi==0:
        return None
    L,H=max(0,S.lo),min(limit,S.hi)
    if state=='A':return Box(L,H),B0
    if state=='B':return B0,Box(L,H)
    assert state in ('X','F')
    val=Box(max(0,L-SCALE),min(SCALE,H))
    return val,val


def early_code_tree(j,word,mid,support_masks,out_count):
    initial_a,initial_b=setup(t_box())
    a,b=[point(x) for x in initial_a],[point(x) for x in initial_b]
    t=sequence(j,(word,mid,()))
    U=[(1,)]*4+support_masks
    expand=[1]*(len(U)+1)
    for k in range(len(U)-1,-1,-1):
        expand[k]=expand[k+1]*len(U[k])
    def path(k,bits):
        if k==len(U):
            yield (word,mid,tuple(bits)),tuple(a),tuple(b)
            return
        state=t[k]
        internal=cross(a,b,k+8)
        for C in U[k]:
            # Only the designated early positive rows have strict support.
            selected=fill_interval(C-internal,state,k<len(word) and state!='Z')
            if selected is None:
                out_count['early_interval']+=expand[k+1]
                continue
            A,B=selected
            a.append(A);b.append(B)
            if k>=4:bits.append(C)
            yield from path(k+1,bits)
            if k>=4:bits.pop()
            a.pop();b.pop()
    return path(0,[])


def census(j):
    tally=Counter()
    late=[]
    early_cert=[]
    total=2**(j-19)
    for early,mid in initial_shapes(j):
        fixed=masks(j,early)
        tally['all']+=total
        tally['positive_support_mask']+=total-2**sum(len(d)==2 for d in fixed)
        for shape,A,B in early_code_tree(j,early,mid,fixed,tally):
            positions=slots(j,shape)
            if positions:
                E=exact_system(j,shape)
                coeff=[f for _,f in E if all(len(a)<2 for a in f)]
                n=len(positions)
                M=[tuple(f.get((i,),ZERO) for i in range(n)) for f in coeff]
                rhs=[-f.get((),ZERO) for f in coeff]
                if len(coeff)==n:
                    witness=linear_test(positions,M,rhs)
                    if witness:
                        tally['norm']+=1
                        early_cert.append((shape,'norm',witness))
                        continue
                quick=local_linear_rejection(M,rhs)
                if quick:
                    tally[quick[0]]+=1
                    early_cert.append((shape,*quick))
                    continue
            tally['retained']+=1
            late.append({'shape':shape,'A':[(v.lo,v.hi) for v in A],'B':[(v.lo,v.hi) for v in B]})
    assert tally['all']==sum(v for k,v in tally.items() if k!='all')
    return {'degree':j,'counts':dict(tally),'late':late,'linear_certificates':early_cert}


def late_codes(j,shape,A,B,ledger):
    a,b=[Box(*w) for w in A],[Box(*w) for w in B]
    ap,bp=required_positive(shape[0])
    suffix=[(1,) if proof_exists_positive(ap,bp,k) else (0,1) for k in range(j-7,j)]
    ledger['support_word_count']+=2**sum(len(z)==2 for z in suffix)
    future=[1]*8
    for i in range(6,-1,-1):future[i]=future[i+1]*len(suffix[i])
    def branch(i,word):
        if i==7:
            yield tuple(word)
            return
        k=j-7+i
        side='Z' if k==j-1 else 'B' if k in (j-2,j-5) else 'A' if k in (j-3,j-4) else 'F'
        internal=cross(a,b,k)
        for C in suffix[i]:
            selected=fill_interval(C-internal,side)
            if selected is None:
                ledger['tail_interval_rejected']+=future[i+1]
                continue
            x,y=selected
            a.append(x)
            b.append(y)
            word.append(C)
            yield from branch(i+1,word)
            word.pop();a.pop();b.pop()
    return list(branch(0,[]))


if __name__=='__main__':
    if not __debug__:raise SystemExit('requires assertions')
    D=[]
    for j in map(int,sys.argv[1:]):
        r=census(j)
        D.append(r)
        tail=Counter()
        for item in r['late']:
            words=late_codes(j,item['shape'],item['A'],item['B'],tail)
            if words:tail['tail_live_shapes']+=1
            tail['tail_live_codes']+=len(words)
        print(j,r['counts'],dict(tail),flush=True)
    Path('prefix'+'.'.join(sys.argv[1:])+'.json').write_text(json.dumps(D,separators=(',',':'))+'\n')
