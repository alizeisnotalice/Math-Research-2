#!/usr/bin/env python3
"""Exact finite tensor orbit/JD grouping ledger; no LP or MATLAB launch."""
from pathlib import Path
from fractions import Fraction as F
import argparse,hashlib,json,os,sys,time

def main():
    assert os.path.ismount('/Volumes/data')
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--grid',required=True);ap.add_argument('--output',required=True);args=ap.parse_args()
    started=time.monotonic();p=Path(args.grid);grid=json.loads(p.read_text());e=list(map(F,grid['axis_edges']));m=len(e)-1
    assert m%2==1 and all(e[i]<e[i+1]for i in range(m)) and all(e[k]==-e[m-k]for k in range(m+1))
    mid=m//2;index={};J=0
    for i in range(mid,m):
        for j in range(mid,i+1):
            J+=1
            for ii in (i,m-1-i):
                for jj in (j,m-1-j):
                    for a,b in [(ii,jj),(jj,ii)]:
                        assert (a,b)not in index or index[a,b]==J
                        index[a,b]=J
    assert len(index)==m*m
    rows=[]
    for axis in (0,1):
        for k in range(m+1):
            for t in range(m):
                plus=index.get((k,t)if axis==0 else(t,k),0)
                minus=index.get((k-1,t)if axis==0 else(t,k-1),0)
                if plus==minus:pair=None;sign=0
                else:pair=tuple(sorted((plus,minus)));sign=1 if plus==pair[1]else-1
                physical={}
                if plus:physical[plus]=physical.get(plus,0)+1
                if minus:physical[minus]=physical.get(minus,0)-1
                physical={q:a for q,a in physical.items()if a}
                canonical={}if pair is None else{pair[1]:1,**({pair[0]:-1}if pair[0]else{})}
                assert physical=={q:sign*a for q,a in canonical.items()}
                rows.append((pair,sign))
    pairs=sorted({pair for pair,sign in rows if sign});lookup={pair:i+1 for i,pair in enumerate(pairs)};G=len(pairs)
    group=[lookup[pair]if sign else 0 for pair,sign in rows];signs=[sign for pair,sign in rows];counts=[group.count(k+1)for k in range(G)]
    out=dict(status='EXACT_INTEGER_PHYSICAL_JD_EQUALS_SIGNED_GROUP_MAP',scope='Exact finite geometry/group-map construction only; no LP feasible-set comparison, candidate or source residual computation.',python=sys.version,script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),grid=str(p.resolve()),grid_sha256=hashlib.sha256(p.read_bytes()).hexdigest(),axis_intervals=m,physical_cells=m*m,D4_orbits=J,physical_interfaces=len(rows),zero_interface_count=group.count(0),canonical_jump_groups=G,group_multiplicity_range=[min(counts),max(counts)],subdiv2_target_subboxes=4*J,subdiv2_grouped_variables=J+1+2*G+8*J,subdiv2_old_variables=J+1+2*len(rows)+8*J,group_pairs=pairs,physical_group_1based=group,physical_sign=signs,elapsed_seconds=time.monotonic()-started)
    dest=Path(args.output);assert str(dest.resolve()).startswith('/Volumes/data/');dest.parent.mkdir(parents=True,exist_ok=True);dest.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items()if k not in ['group_pairs','physical_group_1based','physical_sign']},indent=2))
if __name__=='__main__':main()
