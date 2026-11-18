"""Exact prefix arithmetic and certified finite-field eliminants for degrees 25/26."""
from argparse import ArgumentParser
from collections import Counter
from functools import lru_cache
from hashlib import sha256
from math import isqrt
from itertools import product
from pathlib import Path
import json,sys,time

import sympy as sp
from sympy.polys import distributedmodules as dm
from sympy.polys.orderings import lex
from sympy.polys.rings import ring

sys.dont_write_bytecode=True
from verify_degrees_23_24 import T,ZERO,setup,positive_low_rows,proof_exists_positive,linear_test,local_linear_rejection,INTERVAL,H,evaluate
from verify_degree_19 import cross
from verify_degrees_17_18 import common_values

CODES='ABXZ'
ALLOWED={'A':'A','B':'B','X':'Z','Z':'F'}


def generate(j):
    for x in product('ABX','ABX','ABX',CODES,CODES):
        early=''.join(x)
        for middle in ('AB' if j%2==0 else (None,)):
            for bits in product((0,1),repeat=j-19):
                yield early,middle,bits


def seq(j,shape):
    early,middle,_=shape
    p=early+(middle or '')+''.join(ALLOWED[a] for a in early[::-1])
    assert len(p)==j-15
    return p


def counts(j,shape):
    return tuple(i for i,a in enumerate(seq(j,shape),8) if a in ('X','F'))


def basic_rejection(shape):
    early,mid,bits=shape
    a,b=positive_low_rows()
    for i,p in enumerate(early,8):
        if p in ('A','X'):a.add(i)
        if p in ('B','X'):b.add(i)
    return any(not c and proof_exists_positive(a,b,k) for k,c in enumerate(bits,12))


def model(t,j,shape,vars):
    z=iter(vars)
    res=[]
    A,B=setup(t)
    for k,(side,c) in enumerate(zip(seq(j,shape),[1]*4+list(shape[2]),strict=True),8):
        r=c-cross(A,B,k)
        if side=='A':a,b=r,0
        elif side=='B':a,b=0,r
        elif side in ('X','F'):
            a=next(z)
            b=r-a
        else:
            a=b=0
            res.append((k,r))
        A.append(a);B.append(b)
    assert not list(z)
    return res,A,B


def early_matrix(j,shape):
    slots=counts(j,shape)
    n=len(slots)
    zero=[0]*n
    base=model(T,j,shape,zero)[0]
    assert len(base)==n
    rows=[k for k,_ in base]
    base=[w for _,w in base]
    def f(x):
        return [w for _,w in model(T,j,shape,x)[0]]
    cols=[]
    quads={}
    for i in range(n):
        one=zero[:]
        one[i]=1
        vals=f(one)
        cols.append([a-b for a,b in zip(vals,base,strict=True)])
        one[i]=2
        vals2=f(one)
        # coefficient of x_i^2 multiplied by 2.
        q=[a-2*b+c for a,b,c in zip(vals2,vals,base,strict=True)]
        if any(x!=ZERO for x in q):
            assert all(all(v%2==0 for v in a.coefficients) for a in q)
            q=[type(a)(tuple(v//2 for v in a.coefficients)) for a in q]
            cols[-1]=[x-y for x,y in zip(cols[-1],q,strict=True)]
            quads[i,i]=q
        for k in range(i):
            two=zero[:]
            two[i]=two[k]=1
            one2=zero[:]
            one2[k]=1
            q=[a-b-c+d for a,b,c,d in zip(f(two),vals,f(one2),base,strict=True)]
            if any(x!=ZERO for x in q):
                quads[k,i]=q
    return rows,slots,base,list(zip(*cols,strict=True)),quads


def scout(j):
    c=Counter();live=[];pure=[]
    for shape in generate(j):
        c['all']+=1
        if basic_rejection(shape):
            c['zeros']+=1
            continue
        n=len(counts(j,shape))
        if not n:
            _,A,B=model(INTERVAL,j,shape,[])
            if any(hasattr(q,'hi') and q.hi<0 for q in A+B):c['negative']+=1
            else:c['pure']+=1;pure.append(shape)
            continue
        r,slots,b,M,Q=early_matrix(j,shape)
        linear=[k for k in range(n) if not any(q[k]!=ZERO for q in Q.values())]
        if Q:
            lin=local_linear_rejection([M[k] for k in linear],[-b[k] for k in linear])
            if lin:c[lin[0]]+=1;continue
            counts_q=(len(Q),len(linear),n)
            c['nonlinear']+=1
            live.append((shape,'quadratic',slots,counts_q))
        else:
            rhs=[-v for v in b]
            cert=linear_test(slots,M,rhs)
            if cert:c['norm']+=1;continue
            cert=local_linear_rejection(M,rhs)
            if cert:c[cert[0]]+=1;continue
            c['linear_other']+=1;live.append((shape,'linear',slots,()))
    return {'j':j,'census':dict(c),'live':live,'pure':pure}


MODS=((1279,168,(875,178,686,823,163,1111,1)),(163,-6,(120,112,161,11,78,6,1)))
BITS=tuple(product((0,1),repeat=7))


def check_modulus(prime,t,expected):
    assert all(prime%k for k in range(2,isqrt(prime)+1))
    assert evaluate(H,t)%prime==0
    a,b=setup(t)
    sources=(a,a,a,b,b,a)
    weight=[0]*6+[1]
    for s in range(5,-1,-1):
        weight[s]=-sum(weight[r]*sources[s][r-s] for r in range(s+1,7))%prime
    assert tuple(weight)==expected
    for there,positions in ((b,(0,1,3,4)),(a,(0,1,2,5))):
        for s in positions:
            assert sum(weight[r]*there[r-s] for r in range(s,7))%prime==0
    if prime==1279:
        hashed=Counter(sum(x*y for x,y in zip(weight,d,strict=True))%prime for d in BITS)
        assert len(hashed)==127 and hashed.pop(0)==2 and set(hashed.values())=={1}
        assert sum(weight)%prime==0


@lru_cache
def make_ring(prime,n):
    return ring([f'u{i}' for i in range(n)]+['z'],sp.GF(prime))


def truncated_product(A,B,k):
    return sum(A[i]*B[k-i] for i in range(max(0,k-len(B)+1),min(k,len(A)-1)+1))


def digits_by_support(j,shape):
    a,b=positive_low_rows()
    for k,state in enumerate(shape[0],8):
        if state in ('A','X'):
            a.add(k)
        if state in ('B','X'):
            b.add(k)
    return [word for word in BITS if all(d or not proof_exists_positive(a,b,i)
                                           for i,d in enumerate(word,start=j-7))]


def pack(poly,prime):
    return [[list(exponents),int(coeff)%prime] for exponents,coeff in sorted(poly.items())]


def unpack(poly,R):
    return R.from_dict({tuple(e):R.domain.convert(c) for e,c in poly})


def field_certificate(j,shape,mod,witness=None):
    prime,t,weight=mod
    R,*vars=make_ring(prime,len(counts(j,shape)))
    z=vars[-1]
    all_eq,A,B=model(R(t),j,shape,vars[:-1])
    equations=[R(f) for _,f in all_eq]+[z-sum(w*truncated_product(A,B,k)
                                              for k,w in enumerate(weight,j-7))]
    certs=[]
    if witness is None:
        symbols=tuple(u.as_expr() for u in vars)
        elements=[dm.sdm_from_vector([g.as_expr()],lex,R.domain,gens=symbols) for g in equations]
        basis,transition=dm.sdm_groebner(elements,dm.sdm_nf_mora,lex,R.domain,extended=True)
        for poly,multipliers in zip(basis,transition,strict=True):
            if any(any(m[1:-1]) for m,c in poly):
                continue
            g=R.from_expr(dm.sdm_to_vector(poly,symbols,R.domain,n=1)[0])
            m=[R.from_expr(x) for x in dm.sdm_to_vector(multipliers,symbols,R.domain,n=len(equations))]
            certs.append((g,m))
    else:
        certs=[(unpack(proof['only_z'],R),[unpack(a,R) for a in proof['multipliers']]) for proof in witness]
    results=[]
    univariate=[]
    for f,m in certs:
        assert f and all(not any(e[:-1]) for e in f)
        assert len(m)==len(equations)
        # This is the certificate. The checker does not assume that its
        # generators form a Groebner basis or trust the solver status.
        assert f==sum(a*b for a,b in zip(m,equations,strict=True))
        results.append({'only_z':pack(f,prime),'multipliers':[pack(g,prime) for g in m]})
        univariate.append(tuple((e[-1],int(c)%prime) for e,c in f.items()))
    return results,univariate


def restrict(words,mod,univariate):
    prime,_,weight=mod
    live=[]
    for bits in words:
        z=sum(x*y for x,y in zip(weight,bits,strict=True))%prime
        if all(sum(c*pow(z,k,prime) for k,c in f)%prime==0 for f in univariate):
            live.append(bits)
    return live


def shape_lists(report):
    return [shape for shape,*_ in report['live']]+report['pure']


def canonical_shift(shape):
    a,b,c=shape
    return a,b,tuple(c)


def project(early,witness=None):
    degree=early['j']
    census=Counter()
    all_rows=[]
    there={canonical_shift(row['early_shape']):row for row in witness} if witness is not None else None
    for shape in shape_lists(early):
        words=digits_by_support(degree,shape)
        known=there[shape] if there is not None else None
        proofs={}
        count_start=len(words)
        census['binary_words_start']+=count_start
        for p,t,w in MODS:
            if not words:
                break
            key=str(p)
            orig=known['proofs'][key] if known is not None else None
            certs,univariate=field_certificate(degree,shape,(p,t,w),orig)
            census[f'mod{p}_shape_checks']+=1
            census[f'mod{p}_identities']+=len(certs)
            before=len(words)
            words=restrict(words,(p,t,w),univariate)
            census[f'mod{p}_rejected_words']+=before-len(words)
            proofs[key]=certs
        assert not words,(degree,shape,words)
        all_rows.append({'early_shape':shape,'late_words_from_support':count_start,'proofs':proofs})
    assert len(all_rows)==len(shape_lists(early))
    if there is not None:
        assert len(there)==len(all_rows)
    return dict(census),all_rows


def main():
    if not __debug__:raise SystemExit('run with assertions enabled')
    cli=ArgumentParser(description=__doc__)
    cli.add_argument('--out',type=Path,required=True)
    cli.add_argument('--witness',type=Path)
    args=cli.parse_args()
    start=time.monotonic()
    expected={
        25:dict(all=27648, zeros=26999, negative=54, norm=257, zero_linear_row=146,
                pure=1, local_nonintegral_norm=57, modular_linear_row=4, linear_other=3, nonlinear=127),
        26:dict(all=110592, zeros=109054, negative=162, norm=579, zero_linear_row=422,
                pure=2, linear_other=5, local_nonintegral_norm=114, nonlinear=254),
    }
    orig=None if args.witness is None else {x['degree']:x['certificates']
         for x in json.loads(args.witness.read_bytes())['degrees']}
    for p,t,weights in MODS:check_modulus(p,t,weights)
    degree_reports=[]
    for j in (25,26):
        early=scout(j)
        assert early['census']==expected[j]
        stats,certs=project(early,None if orig is None else orig[j])
        degree_reports.append({'degree':j,'early_counts':expected[j],'modular_counts':stats,
                               'certificates':certs})
        print(j,expected[j],stats,flush=True)
    paths=[Path(__file__).with_name(x) for x in ('verify_degrees_17_18.py','verify_degree_19.py',
                 'verify_degrees_21_22.py','verify_degrees_23_24.py',Path(__file__).name)]
    data={'status':'verified','mode':'generate' if orig is None else 'identities_only',
          'degrees':degree_reports,'unknown':0,'python':sys.version,'sympy':sp.__version__,
          'source_sha256':{p.name:sha256(p.read_bytes()).hexdigest() for p in paths},
          'seconds':time.monotonic()-start}
    if args.witness is not None:data['input_receipt_sha256']=sha256(args.witness.read_bytes()).hexdigest()
    with args.out.open('x') as f:f.write(json.dumps(data,separators=(',',':'),sort_keys=True)+'\n')


if __name__=='__main__':main()
