#!/usr/bin/env python3
"""Freeze 12 deterministic structural strata x six cases; no RNG involved."""
import argparse
import hashlib
import json
import os
from pathlib import Path

NAMES = ['near-compact-square-grid','long-single-row','two-row-strip',
         'chain-of-small-clusters','diagonal-chain','dyadic-square-shells',
         'nested-axis-shells','unequal-density-cluster-bridge',
         'square-frame-with-interior','long-rectangular-grid',
         'heavy-interior-light-grid','corner-clusters-with-connectors']


def unique(seq):
    return list(dict.fromkeys(seq))


def design(stratum,case):
    n = (16,17,18,19,20,20)[case]
    variant = case+1
    if stratum == 0:
        p=[(4*x,4*y) for y in range(4) for x in range(5)][:n]
        # Total mass close to d^2: exact variant straddles compactness.
        span=max(max(t[a] for t in p)-min(t[a] for t in p) for a in (0,1))
        target=span*span+(case-2)*n
        m=[target//n+(i<target%n) for i in range(n)]
    elif stratum == 1:
        p=[(4*i,0) for i in range(n)]
        m=[25+((i+variant)%3)*4 for i in range(n)]
        m[-1]=1 if case%2 else 16
    elif stratum == 2:
        p=[(4*(i//2),4*(i%2)) for i in range(n)]
        m=[16+9*((i+case)%3==0) for i in range(n)]
        m[-1]=1+case
    elif stratum == 3:
        p=[(14*(i//4)+2*(i%2),2*((i%4)//2)) for i in range(n)]
        m=[9+7*((i//4+case)%2) for i in range(n)]
    elif stratum == 4:
        p=[(3*i,3*i+(i%2)) for i in range(n)]
        m=[25+11*((i+case)%4==0) for i in range(n)]
    elif stratum == 5:
        p=[]
        for r in (2,4,8,16,32):
            p.extend([(r,r),(-r,r),(-r,-r),(r,-r)])
        p=p[:n]
        m=[max(1,max(abs(x),abs(y))**2//(3+case%3)) for x,y in p]
    elif stratum == 6:
        p=[]
        for r in (1,2,4,8,16):
            p.extend([(r,0),(0,r),(-r,0),(0,-r)])
        p=p[:n]
        m=[max(1,(abs(x)+abs(y))**2//(2+case%3)) for x,y in p]
    elif stratum == 7:
        left=[(-24+2*x,2*y) for y in (-1,0,1) for x in (0,1)]
        right=[(24+2*x,2*y) for y in (-1,0,1) for x in (0,1)]
        bridges=[(-16+32*i//(n-13),(-1)**i) for i in range(n-12)]
        p=left+right+bridges
        m=[64]*6+[9+case]*6+[1+case%3]*(n-12)
    elif stratum == 8:
        p=[(-8,-8),(0,-8),(8,-8),(8,0),(8,8),(0,8),(-8,8),(-8,0),
           (-4,-8),(4,8),(-8,4),(8,-4)]
        interior=[(x,y) for y in (-4,0,4) for x in (-4,0,4)]
        p+=interior[:n-12]
        m=[25+case]*12+[1+case%3]*(n-12)
    elif stratum == 9:
        p=[(8*(i//4),2*(i%4)) for i in range(n)]
        m=[9+7*((i+case)%5==0) for i in range(n)]
    elif stratum == 10:
        ring=[(x,y) for y in (-8,-4,0,4,8) for x in (-8,-4,0,4,8) if (x,y)!=(0,0)]
        ring.sort(key=lambda p:(max(abs(p[0]),abs(p[1])),p[1],p[0]))
        p=[(0,0)]+ring[:n-1]
        m=[128+32*case]+[1+(i+case)%3 for i in range(n-1)]
    else:
        p=[]
        for a,b in ((-16,-16),(16,-16),(16,16),(-16,16)):
            p.extend([(a,b),(a+2,b),(a,b+2),(a+2,b+2)])
        p+=[(-8,0),(0,0),(8,0),(0,8)][:n-16]
        m=[16+9*((i//4+case)%2) for i in range(16)]+[1+case]*(n-16)
    assert len(p)==n and len(unique(p))==n and min(m)>0
    # The second N20 case changes geometry as well as masses.
    if case==5:
        p=[(2*x+(i%3==0),2*y) for i,(x,y) in enumerate(p)]
        m=[4*x+(i%4) for i,x in enumerate(m)]
        assert len(unique(p))==n
    return {'id':f's{stratum:02}-c{case:02}-n{n}',
            'coordinate_unit':'1/4','mass_unit':'1/16',
            'points_quarters':p,'mass_sixteenths':m,
            'design':'Deterministic structural formula; no random seed or optimization.',
            'parameters':{'stratum':stratum,'case':case,'n':n,'variant':variant}}


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True)
    args=ap.parse_args();assert os.path.ismount('/Volumes/data')
    strata=[{'name':name,'inputs':[design(s,c) for c in range(6)]}
            for s,name in enumerate(NAMES)]
    data={'schema':'structural-source-frontier-input-v1','status':'INPUTS_ONLY_NOT_EXECUTION',
          'threshold':'1','dimensions':2,'geometry':'axis-parallel centered common-radius squares',
          'source_counts':[16,17,18,19,20,20],'strata':strata,
          'generator_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
          'purpose':'High-rank/noncompact cut and compact-block deletion tests; no uniform conclusion.'}
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(data,indent=2)+'\n')
    print(json.dumps({'cases':72,'sha256':hashlib.sha256(args.out.read_bytes()).hexdigest()}))


if __name__=='__main__':main()
