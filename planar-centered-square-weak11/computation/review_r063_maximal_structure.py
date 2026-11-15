#!/usr/bin/env python3
"""Independent exact audit of R63 STRUCTURE, not its 3.5 conjecture.

Halo membership uses sorted actual capture radii, never bbox receiver-union
code. Finite point tests support the separately reviewed localization proof;
they are not an area integration or a general deletion-cost certificate.
"""
from collections import Counter
from fractions import Fraction as F
import hashlib
from itertools import product
import json
import os
from pathlib import Path
import sys
import time

HERE=Path(__file__).resolve().parent
SOURCE=Path('/Volumes/data/CodexProjects/T2/.local/workspaces/pcw-astra-covering/planar-centered-square-weak11/proofs/r063-maximal-compact-block-deletion.md')
R61=HERE/'results/round061-source-union-independent-manifest.json'
R61_SHA='879c71cad764f9919424843ff06e103dde4573a9436f664550233a251c690e9e'


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def halo(points,masses,receiver):
    """Check every possible true capture set by increasing L-infinity radius."""
    by_radius=Counter()
    for (x,y),m in zip(points,masses):
        by_radius[max(abs(x-receiver[0]),abs(y-receiver[1]))]+=m
    captured=F(0)
    for r,m in sorted(by_radius.items()):
        captured+=m
        if captured>4*r*r:return True
    return False


def run_case(points,masses,membership):
    n=len(points);full=(1<<n)-1;rows={};L=F(0);best=[]
    for mask in range(1,full+1):
        ids=[i for i in range(n)if mask&(1<<i)]
        a,b=min(points[i][0]for i in ids),max(points[i][0]for i in ids)
        c,d=min(points[i][1]for i in ids),max(points[i][1]for i in ids)
        span=max(b-a,d-c);mass=sum((masses[i]for i in ids),F(0))
        rows[mask]=(mass,span,a,b,c,d)
        if span*span<mass:
            if mass>L:L,best=mass,[mask]
            elif mass==L:best.append(mask)
    assert L>0 and best
    stats=Counter(cases=1,nonempty_subsets=full,maximum_blocks=len(best))
    # Every captured subset of any square side s>=sqrt(L) has
    # mass<=max(L,d(T)^2)<=s^2. This tests ALL subset threshold regimes.
    for mass,span,*_ in rows.values():
        assert mass<=max(L,span*span)
        stats['square_cap_subset_checks']+=1
    for S in best:
        mass,span,a,b,c,d=rows[S]
        assert mass==L and span*span<L
        outside=[i for i in range(n)if not S&(1<<i)]
        for i in outside:
            px,py=points[i]
            assert not(a<=px<=b and c<=py<=d)
            # Membership of the closed collar J is equivalent to q²<=L.
            q=max(b-px,px-a,d-py,py-c)
            assert q>=0 and q*q>L
            stats['closed_collar_point_checks']+=1
        for T in range(1,full+1):
            if T&S:continue
            tm=rows[T][0];um,us,*_=rows[T|S]
            assert um==L+tm and tm<=us*us-L
            stats['exterior_subset_checks']+=1
        local=[i for i,p in enumerate(points)if max(F(0),a-p[0],p[0]-b,c-p[1],p[1]-d)**2<=L]
        lm=sum((masses[i]for i in local),F(0))
        assert lm<=9*L
        # Restriction preserves L exactly, including all tied maximum S.
        lmask=sum(1<<i for i in local)
        assert S&~lmask==0
        assert max((mm for t,(mm,ss,*_)in rows.items()if t&~lmask==0 and ss*ss<mm),default=F(0))==L
        stats['localized_mass_checks']+=1
        if membership:
            pp=[points[i]for i in local];mm=[masses[i]for i in local]
            residual=[i for i in range(n)if not S&(1<<i)]
            lr=[i for i in local if not S&(1<<i)]
            for receiver in product((F(i,2)for i in range(-6,19)),(F(j,2)for j in range(-6,11))):
                A=halo(points,masses,receiver)
                B=halo([points[i]for i in residual],[masses[i]for i in residual],receiver)
                C=halo(pp,mm,receiver)
                D=halo([points[i]for i in lr],[masses[i]for i in lr],receiver)
                q=max(F(0),a-receiver[0],receiver[0]-b,c-receiver[1],receiver[1]-d)
                insideI=4*q*q<L
                if insideI:assert A==C and B==D
                assert (A and not B)==(C and not D)
                if A and not B:assert insideI
                stats['direct_radius_receiver_checks']+=1
    return stats


def main():
    assert os.path.ismount('/Volumes/data')
    started=time.perf_counter();source_sha=sha(SOURCE)
    assert sha(R61)==R61_SHA
    total=Counter();strata=[]
    grids=[('square_four',[(F(x),F(y))for x in (0,2)for y in (0,2)]),
           ('two_row_strip',[(F(x),F(y))for x in (0,2,4,6)for y in (0,1)]),
           ('square_nine',[(F(x),F(y))for x in (-2,0,2)for y in (-2,0,2)])]
    for name,points in grids:
        subtotal=Counter()
        for index,masses in enumerate(product((F(1),F(2)),repeat=len(points))):
            counts=run_case(points,masses,membership=index in (0,1,3,(1<<len(points))-1))
            subtotal.update(counts)
        total.update(subtotal);strata.append(dict(name=name,points=[[str(x),str(y)]for x,y in points],mass_menu=['1','2'],statistics=dict(subtotal)))
    # Boundary/equality examples and arbitrary tiny positive masses are explicit.
    extras=[('singleton',[(F(0),F(0))],[F(1)]),
            ('strict_equal_pair',[(F(0),F(0)),(F(1),F(0))],[F(1,2),F(1,2)]),
            ('tiny_bridge',[(F(0),F(0)),(F(1),F(0)),(F(2),F(0)),(F(3),F(0))],[F(1),F(1,1000),F(1,1000),F(1)]),
            ('diagonal_wide',[(F(0),F(0)),(F(1),F(1)),(F(2),F(2)),(F(4),F(4))],[F(1),F(1),F(2),F(1)])]
    for name,points,masses in extras:
        counts=run_case(points,masses,True);total.update(counts)
        strata.append(dict(name=name,points=[[str(x),str(y)]for x,y in points],masses=list(map(str,masses)),statistics=dict(counts)))
    assert sha(SOURCE)==source_sha and sha(R61)==R61_SHA
    # The R61 five payload files are frozen as well as its manifest.
    old=json.loads(R61.read_text())
    assert len(old['files'])==5
    for row in old['files']:assert sha(Path(row['source']))==row['sha256']
    out=dict(status='PASS_R63_STRUCTURE_ONLY',python=sys.version,executable=sys.executable,script_sha256=sha(Path(__file__)),
             arithmetic='Exact Fraction; no square roots, floating decisions, randomization, bbox halo-union imports or author code calls.',
             proof_source=str(SOURCE),proof_sha256=source_sha,statistics=dict(total),strata=strata,
             R61_manifest_frozen_sha256=R61_SHA,R61_five_payload_hashes_unchanged=True,
             elapsed_seconds=time.perf_counter()-started,
             scope='All-subset square caps, maximal-block closure, closed collar, exterior constraints, local mass/L invariance; finite receiver checks by direct actual capture radii. No area integration or proof/test conclusion about conjecture(6).')
    target=HERE/'results/round063-structure-independent.json';target.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))


if __name__=='__main__':main()
