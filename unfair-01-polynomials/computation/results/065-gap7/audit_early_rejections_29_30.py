"""Recheck every early norm and modular rejection by polynomial resultants."""
from argparse import ArgumentParser
from collections import Counter
from functools import cache
from hashlib import file_digest
import json
from pathlib import Path

from sympy import Poly, ZZ, symbols
from prefix_29_30 import ZERO, canon, exact_system, slots
from verify_degrees_17_18 import H


t=symbols('t')
h=Poly.from_list(H[::-1],t,domain=ZZ)
zero=Poly(0,t,domain=ZZ)
one=Poly(1,t,domain=ZZ)


@cache
def poly(value):
    return Poly.from_list(value.coefficients[::-1],t,domain=ZZ)


def determinant(rows):
    """Laplace expansion over Z[t]/H, retaining the least supported row."""
    if not rows:
        return one
    if len(rows)==1:
        return rows[0][0]
    i=min(range(len(rows)),key=lambda k:sum(bool(x) for x in rows[k]))
    answer=zero
    remaining=rows[:i]+rows[i+1:]
    for j,v in enumerate(rows[i]):
        if v:
            minor=tuple(row[:j]+row[j+1:] for row in remaining)
            term=(v*determinant(minor)).rem(h)
            answer=answer-term if (i+j)%2 else answer+term
    return answer


@cache
def norm(f):
    return int(h.resultant(f))


def check_row(j,record):
    shape,kind,*witness=record
    shape=canon(shape)
    positions=slots(j,shape)
    n=len(positions)
    affine=[f for k,f in exact_system(j,shape) if all(len(monomial)<2 for monomial in f)]
    M=[tuple(poly(f.get((i,),ZERO)) for i in range(n)) for f in affine]
    b=[-poly(f.get((),ZERO)) for f in affine]
    if kind=='norm':
        k,numerator,denominator=witness[0]
        assert len(M)==n
        c=positions.index(k)
        numerator_matrix=tuple(row[:c]+(rhs,)+row[c+1:] for row,rhs in zip(M,b,strict=True))
        assert norm(determinant(tuple(M)))==denominator and denominator!=0
        assert norm(determinant(numerator_matrix))==numerator and numerator%denominator!=0
    elif kind=='local_nonintegral_norm':
        numerator,denominator=witness
        assert denominator!=0 and numerator%denominator!=0
        assert any(len(xs)==1 and norm(xs[0])==denominator and norm(rhs)==numerator
                   for row,rhs in zip(M,b,strict=True) if len(xs:=[x for x in row if x])==1)
    elif kind=='modular_linear_row':
        residue=witness[0]
        assert residue%163!=0
        assert any(all(q.eval(-6)%163==0 for q in row) and rhs.eval(-6)%163==residue
                   for row,rhs in zip(M,b,strict=True))
    else:
        assert kind=='zero_linear_row'
        assert any(not any(row) and rhs for row,rhs in zip(M,b,strict=True))
    return kind


def sha(path):
    with path.open('rb') as f:
        return file_digest(f,'sha256').hexdigest()


def main():
    if not __debug__:
        raise SystemExit('requires assertions')
    parser=ArgumentParser(description=__doc__)
    parser.add_argument('--prefix',type=Path,required=True)
    parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args()
    data=json.loads(args.prefix.read_text())
    counts={}
    for case in data:
        j=case['degree']
        rows=Counter(check_row(j,record) for record in case['linear_certificates'])
        assert dict(rows)=={kind:n for kind,n in case['counts'].items() if kind in ('norm','modular_linear_row','local_nonintegral_norm','zero_linear_row')}
        counts[str(j)]=dict(rows)
        print(j,dict(rows),flush=True)
    result={'status':'verified','counts':counts,'input_sha256':sha(args.prefix),'source_sha256':sha(Path(__file__))}
    with args.out.open('x') as out:
        out.write(json.dumps(result,sort_keys=True,indent=2)+'\n')


if __name__=='__main__':
    main()
