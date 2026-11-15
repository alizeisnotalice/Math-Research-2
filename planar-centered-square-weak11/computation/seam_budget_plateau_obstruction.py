#!/usr/bin/env python3
"""Bounded exact checks for the density-two plateau seam obstruction.

The paired proof controls every radius and the full halo analytically.
This script verifies rational block geometry, strict witnesses, transport
bounds, and area estimates. It neither searches nor computes the full halo.
"""
from fractions import Fraction as F
from pathlib import Path
import hashlib,json,os,platform,sys


def check(a,K,n):
    h=a/n; s=h/2; L=4*a*K; mass=4*a*L
    nx=2*n; ny=4*K*n
    atoms=[]
    for i in range(nx):
        x=-a+(F(i)+F(1,2))*h
        for j in range(ny):
            y=-L/2+(F(j)+F(1,2))*h
            atoms.append((x,y,2*h*h,j//(4*n)))
    assert len(atoms)==8*K*n*n
    assert sum(p[2] for p in atoms)==mass
    # Every point of every grid cell is at distance <=h/2 from its atom.
    # The following rational radius is strict, so all of R is in the halo.
    singleton_witness_radius=3*h/5
    assert singleton_witness_radius>s
    assert 2*h*h>4*singleton_witness_radius**2
    block_right=[]
    for j in range(K):
        block=[p for p in atoms if p[3]==j]
        m=sum(p[2] for p in block)
        assert m==16*a*a
        rc=2*a
        xmin=min(p[0] for p in block); xmax=max(p[0] for p in block)
        ymin=min(p[1] for p in block); ymax=max(p[1] for p in block)
        assert xmax-xmin==2*a-h and ymax-ymin==4*a-h
        assert 4*rc*rc==m
        assert xmax-xmin<2*rc and ymax-ymin<2*rc
        assert xmin+rc==a+h/2
        assert ymin+rc-(ymax-rc)==h
        block_right.append(xmin+rc)
    # At exterior coordinate distance d=h=2s, the polynomial difference
    # (r+s)(r+s-d)-r^2 has zero r-coefficient and negative constant term.
    d=h
    assert 2*s-d==0
    assert s*(s-d)==-h*h/4<0
    c_lo=a+h/2; c_hi=a+h
    assert all(q>=c_hi-h/2 for q in block_right)
    # Reflection about the true c_h has delta in [h/2,h]. Two cell systems
    # lie in a density <=2 enclosing rectangle of width <=4a+2h.
    old_area_lower=2*a*L
    old_area_upper=(2*a+2*h)*(L+2*h)
    new_area_upper=(4*a+4*h)*(L+2*h)
    excess_upper=new_area_upper-2*old_area_lower
    assert excess_upper==4*h*L+8*a*h+8*h*h
    # A near-boundary strict L1 witness and its threshold-critical limit.
    sigma=min(a/4,h/3); rr=2*a-sigma
    captured=2*(2*a)*(2*rr)
    assert captured>4*rr*rr
    assert (a-sigma)-rr==-a
    assert 16*a*a==4*(2*a)**2
    exact_owner=None
    if K>=2:
        P=2*n; Q=4*n+1; B=F(P*Q,2)
        assert Q<=ny
        owner_mass=2*h*h*P*Q
        assert owner_mass==16*a*a+4*a*h
        assert 4*a*a+a*h==h*h*B
        assert (2*a+h/4)**2>4*a*a+a*h>(2*a)**2
        compared=0
        for pp in range(1,nx+1):
            for qq in range(1,ny+1):
                if (pp-1)**2>=2*pp*qq or (qq-1)**2>=2*pp*qq:
                    continue
                assert qq<=2*pp+1
                # sqrt(A)-pp <= sqrt(B)-P, exact squared comparison.
                A=F(pp*qq,2); dd=P-pp
                if not dd:
                    assert A<=B
                    assert (A==B)==(qq==Q)
                else:
                    C=B-A-dd*dd
                    assert C>0 and C*C>4*dd*dd*A
                compared+=1
        # Sliding exact-owner blocks cover every row, including both ends.
        covered=set()
        for start in range(ny-Q+1): covered.update(range(start,start+Q))
        assert covered==set(range(ny))
        exact_owner={'p':P,'q':Q,'mass_per_owner':str(owner_mass),
                     'source_union_mass':str(mass),'bbox_shapes_compared':compared,
                     'c_formula':'-a+h/2+sqrt(4*a^2+a*h)',
                     'c_lower':str(a+h/2),'c_strict_upper':str(a+3*h/4)}
    return {'a':str(a),'K':K,'n':n,'h':str(h),'L':str(L),'atoms':len(atoms),
            'mass':str(mass),'c_lower':str(c_lo),'c_upper':str(c_hi),
            'active_window':str(h/2),'active_source_union_mass':str(mass),
            'old_area_lower':str(old_area_lower),'old_area_upper':str(old_area_upper),
            'new_area_upper':str(new_area_upper),'true_excess_area_upper':str(excess_upper),
            'excess_over_active_mass_upper':str(excess_upper/mass),
            'singleton_strict_radius':str(singleton_witness_radius),
            'L1_strict_witness_sigma':str(sigma),'L1_strict_witness_radius':str(rr),
            'exact_owner':exact_owner}


def main():
    assert os.path.ismount('/Volumes/data')
    records=[check(F(1),K,n) for K,n in [(1,1),(1,2),(1,4),(1,8),(2,1),(2,2),(2,4)]]
    output={'status':'PASS','scope':'Rational geometry and scalar certificates; no full-halo numerical evaluation.',
            'python':sys.version,'executable':sys.executable,'platform':platform.platform(),
            'lambda':'1','arithmetic':'fractions.Fraction','single_process':True,
            'random_seed':None,'floating_tolerance':None,
            'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'cases':records}
    path=Path(__file__).parent/'results/seam-budget-plateau/exact-check.json'
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(output,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({'status':'PASS','cases':len(records),'max_atoms':max(q['atoms'] for q in records),
                      'output':str(path)},ensure_ascii=False))

if __name__=='__main__': main()
