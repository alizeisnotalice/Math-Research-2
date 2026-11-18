"""Bounded degree-31 check of saturation, tail equations, and both overhangs."""
from argparse import ArgumentParser
from collections import Counter
from hashlib import file_digest
from itertools import product
import json
from pathlib import Path
import sys
import time

from sympy import GF, Matrix, __version__ as sympy_version
from sympy.polys.matrices import DomainMatrix

sys.dont_write_bytecode=True
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'round065-gap7-next'))
from prefix_29_30 import (Box, SCALE, canon, early_code_tree, exact_system, fill_interval,
                          late_codes, masks, required_positive, sequence, slots)
from verify_degree_19 import cross, element
from verify_degrees_25_26 import MODS, T, ZERO, proof_exists_positive, setup
from verify_degrees_29_30 import exact_remaining, individual, modvalue

DEGREE=31
SAMPLES=('AAAAAAAA','AABBBBBB','ABXZABXZ','BAXXAZBX','BBBABABA',
         'XXXXAAAA','XXXXXXXX','XXXZZZZZ','XXXZXZXZ')
TAIL_RULE='FFBAABZ'


def endpoint_matrix():
    results=[]
    for p,t,weights in MODS:
        A,B=setup(t)
        columns=[(B,k) for k in (0,1,3,4)]+[(A,k) for k in (0,1,2,5)]
        raw=[[low[i-k]%p if i>=k else 0 for low,k in columns] for i in range(7)]
        X=DomainMatrix.from_Matrix(Matrix(raw).T).convert_to(GF(p))
        rank=X.rank()
        rows=X.nullspace().to_Matrix()
        assert rank==6 and rows.shape==(1,7)
        divisor=pow(int(rows[0,6]),-1,p)
        vector=tuple(int(x)*divisor%p for x in rows.row(0))
        assert vector==weights
        results.append({'p':p,'t':t,'rank':rank,'weights':weights,'matrix':raw})
    return results


def positive_census(early):
    data=[]
    count=Counter()
    support=masks(DEGREE,early)
    count['all']=1<<len(support)
    count['support']=count['all']-(1<<sum(len(w)==2 for w in support))
    for shape,A,B in early_code_tree(DEGREE,early,None,support,count):
        data.append((shape,A,B))
    count['prefix_kept']=len(data)
    assert count['all']==count['support']+count['early_interval']+count['prefix_kept']
    return data,count


def independent_equations(shape,coordinates):
    """Build the whole polynomial product for each prefix, without cross()."""
    digits=[1]*4+list(shape[2]); free=iter(coordinates)
    A,B=map(list,setup(T))
    eq=[]
    for k,C in enumerate(digits,8):
        inside=sum(a*b for i,a in enumerate(A) for m,b in enumerate(B) if i+m==k)
        rel=C-inside
        side=sequence(DEGREE,shape)[k-8]
        if side=='A':row=(rel,0)
        elif side=='B':row=(0,rel)
        elif side in ('F','X'):
            x=next(free)
            row=(x,rel-x)
        else:
            assert side=='Z'
            row=(0,0)
            eq.append((k,element(rel)))
        A.append(row[0]);B.append(row[1])
    assert not list(free)
    return eq


def evaluate(coefficients,coordinates):
    total=ZERO
    for powers,v in coefficients.items():
        for i in powers:
            v=v*coordinates[i]
        total+=v
    return total


def tail_prefix(A,B,word):
    """Preserve the actual interval prefix for each accepted late word."""
    assert len(A)==len(B)==24 and len(word)==7
    A=list(A); B=list(B)
    for row,C in enumerate(word,24):
        internal=sum(a*B[row-i] for i,a in enumerate(A) if 0<=row-i<len(B))
        pair=fill_interval(C-internal,TAIL_RULE[row-24])
        if pair is None:
            return None
        A.append(pair[0]);B.append(pair[1])
    return A,B


def extension(A,B,side,report):
    """Check degrees 32 through 38 in the chosen 31-leading orientation."""
    one,zero=Box(SCALE,SCALE),Box(0,0)
    A=A+[one if side=='A' else zero]
    B=B+[one if side=='B' else zero]
    compulsory=(1,2,5) if side=='A' else (1,3,4)
    other='B' if side=='A' else 'A'
    seen=set()

    def step(row):
        if row==39:
            return True
        seen.add(row)
        X=cross(A,B,row)
        digits=(1,) if row-31 in compulsory else (0,1)
        for d in digits:
            pair=fill_interval(d-X,other)
            if pair is None:
                continue
            A.append(pair[0]);B.append(pair[1])
            survives=step(row+1)
            A.pop();B.pop()
            if survives:
                return True
        return False

    keep=step(32)
    if not keep:
        report[max(seen)]+=1
    return keep


def run():
    started=time.monotonic()
    matrix=endpoint_matrix()
    totals=Counter()
    logs=[]
    remaining=[]
    exact_witness=None
    for early in SAMPLES:
        paths,c=positive_census(early)
        for shape,A,B in paths:
            coeff=exact_system(DEGREE,shape)
            c['quadratic_shapes']+=bool(any(any(len(m)==2 for m in row) for k,row in coeff))
            for shift in (-2,1):
                coordinates=[(3*i+shift)%7-3 for i in range(len(slots(DEGREE,shape))) ]
                assert independent_equations(shape,coordinates)==[(k,evaluate(row,coordinates)) for k,row in coeff]
                c['exact_assignment_checks']+=1
            independent=[]
            ap,bp=required_positive(early)
            permitted=[(1,) if proof_exists_positive(ap,bp,k) else (0,1) for k in range(24,31)]
            for word in product(*permitted):
                tail=tail_prefix(A,B,word)
                if tail is None:
                    continue
                independent.append(word)
                dead_a=Counter(); dead_b=Counter()
                a=extension(*tail,'A',dead_a)
                b=extension(*tail,'B',dead_b)
                c['A_degree31_overhang_left']+=a
                c['B_degree31_overhang_left']+=b
                c['both_overhangs_impossible']+=not a and not b
                for k,v in dead_a.items():c[f'A_first_empty_by_{k}']+=v
                for k,v in dead_b.items():c[f'B_first_empty_by_{k}']+=v
                for mod in MODS:
                    z=modvalue(word,mod)
                    ok,_=individual(DEGREE,shape,mod,z)
                    if ok:
                        c[f'mod{mod[0]}_rejected']+=1
                        if exact_witness is None and slots(DEGREE,shape):
                            proof,coef=individual(DEGREE,shape,mod,z,track=True)
                            assert proof and coef is not None
                            assert individual(DEGREE,shape,mod,z,witness=coef)[0]
                            exact_witness={'shape':shape,'word':word,'modulus':mod[0],'value':z,'multipliers':coef}
                        break
                else:
                    remaining.append(exact_remaining(DEGREE,shape,[word]))
                    c['mod_unresolved']+=1
            assert independent==late_codes(DEGREE,shape,[(x.lo,x.hi) for x in A],[(x.lo,x.hi) for x in B],c)
            c['mod_start_words']+=len(independent)
        totals.update(c)
        logs.append({'early':early,'counts':dict(c)})
        print(early,dict(c),flush=True)
    assert totals['mod_start_words']==totals['mod1279_rejected']+totals['mod163_rejected']+totals['mod_unresolved']
    assert exact_witness is not None
    return {'status':'bounded_pilot_complete','degree':DEGREE,'scope':SAMPLES,'prefix_census':logs,
            'totals':dict(totals),'tail_matrix':matrix,'one_original_equation_witness':exact_witness,
            'unresolved':remaining,'python':sys.version,'sympy':sympy_version,'seconds':time.monotonic()-started}


def main():
    if not __debug__:
        raise SystemExit('assertions required')
    p=ArgumentParser(description=__doc__)
    p.add_argument('--out',type=Path,required=True)
    a=p.parse_args()
    r=run()
    with Path(__file__).open('rb') as f:
        r['source_sha256']=file_digest(f,'sha256').hexdigest()
    with a.out.open('x') as f:
        f.write(json.dumps(r,sort_keys=True,indent=2)+'\n')
    print(r['totals'],flush=True)


if __name__=='__main__':
    main()
