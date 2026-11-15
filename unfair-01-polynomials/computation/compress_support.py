#!/usr/bin/env python3
"""Exact check of a support-type compression proposed by a floating LP."""
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import random
import numpy as np
from scipy.optimize import linprog
import scipy
import sympy as S


def signature(I,J):
    diagonals={}
    for i,u in enumerate(I):
        for j,v in enumerate(J):
            diagonals.setdefault(u+v,[]).append((i,j))
    return [diagonals[k] for k in sorted(diagonals)]


def compress(I,J):
    r,s=len(I),len(J)
    d=r+s-2
    def coord(side,index):
        row=[0]*d
        if index:
            row[(0 if side==0 else r-1)+index-1]=1
        return row
    def add(a,b): return [x+y for x,y in zip(a,b)]
    def sub(a,b): return [x-y for x,y in zip(a,b)]
    def edge(e): return add(coord(0,e[0]),coord(1,e[1]))
    groups=signature(I,J)
    eq=[sub(edge(e),edge(g[0])) for g in groups for e in g[1:]]
    ge=[sub(coord(side,k),coord(side,k-1)) for side,size in [(0,r),(1,s)] for k in range(1,size)]
    ge += [sub(edge(b[0]),edge(a[0])) for a,b in zip(groups,groups[1:])]
    assert all(sum(c*c for c in row)<=4 for row in eq+ge)
    res=linprog(np.ones(d),A_ub=-np.array(ge),b_ub=-np.ones(len(ge)),
                A_eq=np.array(eq) if eq else None,b_eq=np.zeros(len(eq)) if eq else None,
                bounds=[(0,None)]*d,method='highs')
    if not res.success:
        raise ArithmeticError('LP candidate unavailable, no certificate')
    candidates=[(row,0) for row in eq]+[(row,1) for row in ge if abs(np.dot(row,res.x)-1)<1e-7]
    rows=[]; rhs=[]
    for row,b in candidates:
        if S.Matrix(rows+[row]).rank()>len(rows):
            rows.append(row);rhs.append(b)
        if len(rows)==d: break
    if len(rows)!=d: raise ArithmeticError('Could not recover an exact vertex')
    H=S.Matrix(rows);b=S.Matrix(rhs);x=H.inv()*b
    assert all((S.Matrix([row])*x)[0]==0 for row in eq)
    assert all((S.Matrix([row])*x)[0]>=1 for row in ge)
    D=abs(H.det());y=D*x
    assert all(z.is_Integer and z>=0 for z in y)
    bound=math.isqrt(5**d)
    assert all(z<=bound for z in y)
    II=[0]+list(map(int,y[:r-1]));JJ=[0]+list(map(int,y[r-1:]))
    assert signature(II,JJ)==groups
    for k in range(d):
        Hk=H.copy();Hk[:,k]=b
        assert all(sum(v*v for v in Hk.row(i))<=5 for i in range(d))
        assert abs(Hk.det())==y[k]
    return {'I':I,'J':J,'compressed_I':II,'compressed_J':JJ,'d':d,
            'coordinate_bound':bound,'active_matrix':rows,'active_rhs':rhs,
            'determinant_abs':int(D),'groups':groups,'exact_verified':True}


if __name__=='__main__':
    rng=random.Random(20261115)
    cases=[([0,1,2,4,9],[0,1,3,4,6,11]),([0,2,4,8,18],[0,2,6,8,12,22])]
    cases += [([0]+sorted(rng.sample(range(1,1000),5)),[0]+sorted(rng.sample(range(1,1000),5))) for _ in range(8)]
    print(json.dumps({'seed':20261115,'sympy':S.__version__,'scipy':scipy.__version__,
                      'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                      'cases':[compress(I,J) for I,J in cases]},indent=2))
