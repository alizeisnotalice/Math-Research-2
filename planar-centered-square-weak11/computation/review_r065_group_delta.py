#!/usr/bin/env python3
"""Exact independent canonical-JD factorization and signed-aggregation audit."""
from collections import Counter,defaultdict
from fractions import Fraction as F
import hashlib,json,os
from pathlib import Path
import sys,time

HERE=Path(__file__).resolve().parent
OWNER=Path('/Volumes/data/CodexProjects/T2/.local/workspaces/pcw-sol-kernel-sources/planar-centered-square-weak11')
V2=OWNER/'computation/cauchy_occupation_tensor_jump_lp_v2.m'
V3=OWNER/'computation/cauchy_occupation_tensor_jump_lp_v3.m'
GRID=OWNER/'computation/inputs/r065-tensor-occupation/central-fine-outer-coarse.json'
LEDGER=OWNER/'computation/results/r065-tensor-occupation/exact-group-map.json'
PROOF=OWNER/'proofs/r065-canonical-jump-group-interface.md'

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def pin(p):return dict(path=str(p),sha256=sha(p),bytes=p.stat().st_size)


def main():
    assert os.path.ismount('/Volumes/data')
    start=time.perf_counter()
    assert sha(V3)=='6121279c0386bc614bec6c949a1451016954693c3486dd9563e341f377a8d376'
    assert sha(PROOF)=='36b61d980fc6f3c6e6833fdaa2b1ecc53985362282df2bd8a9842748866bacbc'
    a,b=V2.read_text(),V3.read_text();anchor='function out=lambda_only_diagnostic'
    assert a[a.index(anchor):]==b[b.index(anchor):] # Helpers/diagnostics unchanged.
    tokens=json.loads(GRID.read_text())['axis_edges'];e=list(map(F,tokens));n=len(e)-1;mid=n//2
    assert n==45 and e==[-x for x in reversed(e)]
    # Independent closed-form D4 triangular orbit ID, no author geometry import.
    def oid(i,j):
        aa,bb=sorted((abs(i-mid),abs(j-mid)),reverse=True)
        return aa*(aa+1)//2+bb+1
    J=max(oid(i,j)for i in range(n)for j in range(n));assert J==276
    rows=[];oriented=[]
    for horizontal in (False,True):
        for boundary in range(n+1):
            for transverse in range(n):
                if horizontal:
                    right=oid(transverse,boundary)if boundary<n else 0
                    left=oid(transverse,boundary-1)if boundary else 0
                else:
                    right=oid(boundary,transverse)if boundary<n else 0
                    left=oid(boundary-1,transverse)if boundary else 0
                row=Counter()
                if right:row[right]+=1
                if left:row[left]-=1
                row={k:v for k,v in row.items()if v}
                rows.append(row);oriented.append((right,left))
    keys=sorted({tuple(sorted(pair))for pair in oriented if pair[0]!=pair[1]})
    gid={pair:i+1 for i,pair in enumerate(keys)};signs=[];groups=[]
    for row,(right,left)in zip(rows,oriented):
        if right==left:assert not row;groups.append(0);signs.append(0);continue
        pair=tuple(sorted((right,left)));lo,hi=pair;sgn=1 if right==hi else -1
        GD={hi:1}
        if lo:GD[lo]=-1
        assert row=={k:sgn*v for k,v in GD.items()}
        signs.append(sgn);groups.append(gid[pair])
    ledger=json.loads(LEDGER.read_text())
    assert ledger['grid_sha256']==sha(GRID)
    assert ledger['group_pairs']==[list(k)for k in keys]
    assert ledger['physical_group_1based']==groups and ledger['physical_sign']==signs
    assert len(rows)==4140 and len(keys)==529 and signs.count(0)==0
    counts=Counter(groups);assert min(counts.values())==4 and max(counts.values())==8
    B=J*4;variables=J+1+2*len(keys)+2*B;assert variables==3543
    # Independent exact interval aggregation with mixed physical orientations.
    # Deterministic rational intervals here test algebra, not actual kernel data.
    interval=[(-F(r%7,5),-F(r%7,5)+1+F(r%3,11))for r in range(len(rows))]
    scenarios=[]
    for name,density in [('constant',[F(0)]+[F(1)]*J),('nonmonotone',[F(0)]+[F((j*7)%19,23)for j in range(1,J+1)])]:
        d=[density[hi]-density[lo]for lo,hi in keys]
        totals=[]
        for direction in (0,1):
            lg=[F(0)]*len(keys);ug=[F(0)]*len(keys);psi=[F(0)]*len(keys);physical=F(0)
            for r in range(direction*len(rows)//2,(direction+1)*len(rows)//2):
                if groups[r]==0:continue
                g=groups[r]-1;sgn=signs[r];l,u=interval[r]
                z=l+F(r%3,2)*(u-l)
                assert l<=z<=u
                if sgn>0:lg[g]+=l;ug[g]+=u
                else:lg[g]-=u;ug[g]-=l
                psi[g]+=sgn*z
                physical+=sum(density[k]*v for k,v in rows[r].items())*z
            assert physical==sum(q*v for q,v in zip(d,psi))
            upper=F(0)
            for g in range(len(keys)):
                assert lg[g]<=psi[g]<=ug[g]
                extra=F(g%4,31);dp=max(F(0),d[g])+extra;dn=max(F(0),-d[g])+extra
                assert dp-dn==d[g]
                upper+=dp*ug[g]-dn*lg[g]
            assert physical<=upper
            totals.append(dict(direction='y'if direction else'x',physical_sum=str(physical),group_upper=str(upper)))
        scenarios.append(dict(name=name,zero_actual_groups=sum(q==0 for q in d),totals=totals))
    # Negative orientation must reverse interval endpoints.
    l1,u1,l2,u2=F(1),F(2),F(3),F(5)
    group_range=(l1-u2,u1-l2);assert group_range==(F(-4),F(-1))
    assert u1-l2>u1-u2 # naive sum of upper endpoints would be unsafe.
    out=dict(status='PASS_R65_V3_GROUP_DELTA_INDEPENDENT',python=sys.version,executable=sys.executable,
             script_sha256=sha(Path(__file__)),elapsed_seconds=time.perf_counter()-start,
             physical_rows=len(rows),groups=len(keys),zero_symbolic_rows=0,multiplicity_histogram=dict(sorted(Counter(counts.values()).items())),
             grouped_variables=variables,all4140integer_factorization_rows_matched=True,all_author_saved_maps_matched=True,
             signed_interval_scenarios=scenarios,negative_orientation_counterexample=dict(correct_range=['-4','-1'],wrong_upper='-3'),
             source_pins=[pin(p)for p in (V2,V3,GRID,LEDGER,PROOF)],
             scope='Exact independent JD=SG*GD and finite signed group interval algebra, unchanged helper diagnostics verified. No MATLAB/LP solve, floating LP feasible-set/optimum equivalence, optimized source/fullD/newC2 certificate.')
    target=HERE/'results/round065-v3-group-independent.json';target.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items()if k not in ('source_pins','signed_interval_scenarios')},indent=2))


if __name__=='__main__':main()
