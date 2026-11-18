"""Independently rebuild 3x3 gap-eight branches from convolution."""
from argparse import ArgumentParser
from collections import Counter
from hashlib import sha256
from itertools import product
from pathlib import Path
import json
import sys
import sympy as sp


t,z,z_second=sp.symbols('t z z_second')


def model(side2,side3):
    A,B=[1,t],[1,1-t]
    residuals={}
    for k,policy in enumerate((side2,side3,'A',side3,side2,'both'),start=2):
        X=sp.expand(sum(A[j]*B[k-j] for j in range(1,k)))
        remaining=1-X
        if policy=='both' and k<4:
            split=z_second if k==3 and side2=='both' else z
            A.append(split)
            B.append(remaining-split)
        elif policy=='A':
            A.append(remaining)
            B.append(0)
        elif policy=='B':
            A.append(0)
            B.append(remaining)
        else:
            A.append(0)
            B.append(0)
            residuals[k]=sp.expand(-remaining)
    return A,B,residuals


def main():
    if not __debug__:raise SystemExit('run without -O')
    cmd=ArgumentParser(description=__doc__)
    cmd.add_argument('--source',required=True,type=Path)
    cmd.add_argument('--primary',required=True,type=Path)
    cmd.add_argument('--out',required=True,type=Path)
    args=cmd.parse_args()
    checks={}
    systems={}
    for j,k in product(('A','B','both'),repeat=2):
        A,B,eq=model(j,k)
        systems[j,k]=eq
        # These two reject by a single term C7 less than unity.
        if (j,k) in (('both','both'),('both','A')):
            expected=A[4]*B[3] if k=='both' else A[5]*B[2]
            assert sp.expand(eq[7]+1-expected)==0
            checks[j+' '+k]='single_product'
    source=json.loads(args.primary.read_bytes())['primitive_nonmonic_coefficients_descending']
    polys={name:sp.Poly.from_list(values,t).as_expr() for name,values in source.items()}
    def res(key,d):
        eq=systems[key]
        return sp.expand(sp.resultant(eq[d],eq[7],z))
    J=t**3*(t-1)*(t-2)-(1-t+t*t)
    assert sp.expand(res(('both','B'),6)+t**4*polys['F17'])==0
    assert sp.expand(res(('A','both'),5)-t*polys['F10'])==0
    assert sp.expand(res(('B','both'),5)-t**2*(t*t+2)*(t*t-2*t+2)*J)==0
    checks.update({'both B':'F17','A both':'F10','B both':'J<0'})
    expected={
        ('A','A'):-t*((1-t)**2+1)*((1-t)**4+1),
        ('A','B'):-t*t*polys['F5'],
        ('B','A'):(t-1)*(t*t+1)*(t**4+1),
        ('B','B'):-t*polys['F6'],
    }
    for key,p in expected.items():
        assert systems[key].keys()=={7}
        assert sp.expand(systems[key][7]-p)==0
        checks[' '.join(key)]=str(sp.factor(p))
    assert len(checks)==9
    modp={'F17':37,'F10':13,'F5':7,'F6':5}
    for name,prime in modp.items():
        poly=sp.Poly(polys[name],t)
        content,primitive=poly.primitive()
        assert content==1 and primitive.LC()>1 and primitive.LC()%prime!=0
        assert sp.Poly(primitive,t,modulus=prime).is_irreducible
    data={'status':'verified','source_sha256':sha256(args.source.read_bytes()).hexdigest(),
          'primary_witness_sha256':sha256(args.primary.read_bytes()).hexdigest(),
          'coverage_3x3':checks,'irreducible_mod_prime':modp,'sympy':sp.__version__,'python':sys.version}
    with args.out.open('x') as f:f.write(json.dumps(data,sort_keys=True,indent=2)+'\n')
    print(json.dumps(data,sort_keys=True))


if __name__=='__main__':main()
