#!/usr/bin/env python3
"""Strict per-root summary for the fixed-root r021 production batches."""
import csv, json, math
from decimal import Decimal as D, localcontext
from pathlib import Path

HERE=Path(__file__).resolve().parent
N=121; P=65536; R=64; L=16; T=12
COORDS=[(0,0),(1,0),(1,1),(2,1),(2,2),(3,2),(3,3),(4,3),(4,4),(5,4),(5,5)]
ROOTS=[y*11+x for x,y in COORDS]
assert all(abs(bx-ax)+abs(by-ay)==1
           for (ax,ay),(bx,by) in zip(COORDS,COORDS[1:]))

def summary(root, expected_xy):
    base=f"n11-root{root}-p65536-r64-l16.csv"
    corrected=f"n11-corrected-root{root}-p65536-r64-l16.csv"
    path=HERE/(corrected if (HERE/corrected).exists() else base)
    meta=json.loads(Path(str(path)+'.metadata.json').read_text())
    with path.open() as f: rows=list(csv.DictReader(f))
    assert tuple(int(meta[k]) for k in ('n','N','root','P','R','L','threads')) == (11,121,root,65536,64,16,12)
    assert tuple(meta['root_xy']) == expected_xy and root == expected_xy[1]*11+expected_xy[0]
    assert len(rows)==R and {int(r['rep']) for r in rows}==set(range(R))
    base=int(meta['seed_base'])
    with localcontext() as c:
        c.prec=80
        zs=[]
        for row in rows:
            rep=int(row['rep'])
            assert int(row['root'])==root and int(row['n'])==11 and int(row['N'])==N and int(row['L'])==L
            assert int(row['seed']) == (base+root*1000000007+rep*1009+L*1000003)%(1<<64)
            assert int(row['proposals'])==P*L*(N-2) and int(row['dead'])==0
            assert 0 <= int(row['accept']) <= int(row['proposals'])
            z=D(row['Z']); assert z.is_finite() and z>0; zs.append(z)
        total=sum(zs); mean=total/R
        se=(sum((z-mean)**2 for z in zs)/D(R*(R-1))).sqrt()
        return {'root':root,'xy':[root%11,root//11],'file':path.name,'metadata':meta,
                'mean_h_root':str(mean),'delta_SE':str(se),
                'replicate_weight_ESS':str(total*total/sum(z*z for z in zs)),
                'max_replicate_weight_share':str(max(zs)/total),
                'minimum_particle_parent_ESS':min(float(r['min_parent_ess']) for r in rows),
                'accepted':sum(int(r['accept']) for r in rows),
                'proposals':sum(int(r['proposals']) for r in rows)}

out={'scope':'Finite per-root SMC estimates. Ratios use independent seed batches and delta diagnostics only; no general bound or mixing conclusion.',
     'coordinate_stair':COORDS,
     'runs':[summary(r,xy) for r,xy in zip(ROOTS,COORDS)], 'stair_edges':[]}
byroot={x['root']:x for x in out['runs']}
with localcontext() as c:
    c.prec=80
    for a,b in zip(ROOTS,ROOTS[1:]):
        A,B=byroot[a],byroot[b]
        ma,mb=D(A['mean_h_root']),D(B['mean_h_root']); sa,sb=D(A['delta_SE']),D(B['delta_SE'])
        ratio=mb/ma
        # Independent batch delta approximation for the quotient.
        se=ratio*((sa/ma)**2+(sb/mb)**2).sqrt()
        x,y=A['xy']
        dx,dy=B['xy'][0]-x,B['xy'][1]-y
        assert (dx,dy) in ((1,0),(0,1))
        k=x if dx else y
        # Main candidates indexed by the coordinate changing k -> k+1.
        c_factor=D((k+1)**2)/D(N)
        c_req=ratio*c_factor
        strong_bound=D((11-1-k)**2)/D((k+1)**2)
        # Kept only as a supplementary older stair normalization.
        old_req=ratio*D((x+1)*(y+1))/D(N)
        out['stair_edges'].append({'from_root':a,'to_root':b,'from_xy':[x,y],'to_xy':B['xy'],
                                   'ratio_h_to_over_h_from':str(ratio),'delta_SE':str(se),
                                   'changed_coordinate_k':k,
                                   'main_C_required_for_ratio_le_C_n_squared_over_kplus1_squared':str(c_req),
                                   'main_C_delta_SE':str(se*c_factor),
                                   'strong_bound_nminus1minus_k_over_kplus1_squared':str(strong_bound),
                                   'strong_ratio_over_bound':str(ratio/strong_bound),
                                   'strong_ratio_over_bound_delta_SE':str(se/strong_bound),
                                   'supplementary_old_C_required_xy_product':str(old_req)})
(HERE/'n11-root-summary.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
