from pathlib import Path
import sys
import json
from fractions import Fraction as F
from itertools import combinations
from hashlib import sha256
from flint import arb, ctx
import numpy as np

def sha(path):
    return sha256(path.read_bytes()).hexdigest()

def run(comp, bankdir, out):
    ctx.prec = 192
    sys.path.insert(0, str(comp.resolve()))
    from certify_finite_point_leaf import point_tables, Q, U
    original = json.loads((comp/'finite-point-leaf-input.json').read_bytes())
    cloud = json.loads((bankdir/'new30.json').read_bytes())
    assert sha(bankdir/'new30.json') == "1afef7eabda436e4950853fc4623aa2b3956b625b747bdd7053732fecfe5eb58"
    cells = json.loads((comp/'finite-point-global-input.json').read_bytes())['cells']
    points = cloud['points']
    directions = cloud['directions']
    chosen = [1,2,3,4,5,6,7,8,13,14,25,26,27,28,29]
    small = dict(cloud, points=[points[i] for i in chosen])
    _, ptr, flat = point_tables({**original, "cells":cells})
    wrows = [flat[k:k+8] for k in range(0,len(flat),8)]
    proposals = []
    for j, cell in enumerate(cells):
        ang = float((2*j+1)*arb.pi()/384)
        sn, cs = np.sin(ang), np.cos(ang)
        row = []
        for item in cell:
            ids = item[:3]; labels=[]
            for i in ids:
                a,b=directions[i]
                nx,ny=a*cs+b*sn,-a*sn+b*cs
                labels.append(max(range(len(chosen)),key=lambda k:nx*small['points'][k][0]+ny*small['points'][k][1]))
            row.append(ids+labels)
        proposals.append(row)
    diam, ptr2, updated = point_tables({**small, "cells":proposals})
    assert diam<=small["denominator"]**2 and ptr==ptr2
    assert all(flat[k::8]==updated[k::8] for k in range(7))
    gains=[b-a for a,b in zip(flat[7::8],updated[7::8])]
    assert min(gains)>0
    # A common rational contraction preserves every archived lower score.
    old_scores, new_scores = flat[7::8], updated[7::8]
    assert min(new_scores) > 0
    scale = max(F(a,b) for a,b in zip(old_scores,new_scores))
    assert 0 < scale < F('.999999')
    squashed = {**small, "points":[[scale.numerator*x,scale.numerator*y] for x,y in small["points"]],
                "denominator":scale.denominator*small["denominator"]}
    sqdiam, ptr3, shrunken = point_tables({**squashed, "cells":proposals})
    assert ptr3 == ptr and all(shrunken[k::8]==flat[k::8] for k in range(7))
    assert min(a-b for a,b in zip(shrunken[7::8],old_scores))>=0
    assert F(sqdiam,squashed["denominator"]**2) < F(".999999")**2
    # Divide by a diametral vector as a complex number to make the exact unit input.
    qpts=small["points"]
    farthest=max(((x-u)**2+(y-v)**2,x-u,y-v) for (x,y),(u,v) in combinations(qpts,2))
    sq, ax, ay=farthest
    unit={"points":[[ax*x+ay*y,ax*y-ay*x] for x,y in qpts],"denominator":sq}
    assert max((x-u)**2+(y-v)**2 for (x,y),(u,v) in combinations(unit["points"],2))==sq**2
    orbits=[]
    for i, (x,y) in enumerate(points):
        if y<0: continue
        orbits.append([i] if y==0 else sorted([i,points.index([x,-y])]))
    assert len(orbits)==16 and sorted(i for v in orbits for i in v)==list(range(30))
    indices=np.array([row[:3] for row in wrows])
    ww=np.array([row[3:6] for row in wrows],dtype=float)/np.array([row[6] for row in wrows])[:,None]
    target=np.array([row[7] for row in wrows])/Q
    v=np.array(directions,dtype=float)[indices]/U
    pp=np.asarray(points)/cloud['denominator']
    cr=np.repeat(np.arange(len(cells)),np.diff(ptr))
    guide = []
    for numerator in (0,1,2):
        t=(cr+numerator/2)*2*np.pi/384
        cs,sn=np.cos(t)[:,None],np.sin(t)[:,None]
        nx,ny=v[:,:,0]*cs+v[:,:,1]*sn,-v[:,:,0]*sn+v[:,:,1]*cs
        projections=ww[:,:,None]*(nx[:,:,None]*pp[:,0]+ny[:,:,None]*pp[:,1])
        guide.append(np.stack([np.max(projections[:,:,orb],axis=2) for orb in orbits],axis=-1))
    guide=np.array(guide)
    def row_endpoints(mask):
        subset=[i for i in range(len(orbits)) if mask & (1<<i)]
        return np.take(guide,subset,axis=3).max(axis=3).sum(axis=2) - target
    blockers=[]
    tests={}
    scales=2**48
    def exact_table(part,k):
        key=(part,k)
        if key in tests: return tests[key]
        j=int(cr[k])
        sn,cs=(arb(2*j+part)*arb.pi()/384).sin_cos()
        entries=[]
        row=wrows[k]
        for n,w in zip(row[:3],row[3:6]):
            a,b=directions[n]
            nx,ny=a*cs+b*sn,-a*sn+b*cs
            factor=arb(w)*scales/(row[6]*U*cloud['denominator'])
            upper=[int(((nx*x+ny*y)*factor).upper().ceil().unique_fmpz()) for x,y in points]
            entries.append([max(upper[i] for i in orb) for orb in orbits])
        tests[key]=(entries,row[7]*(scales//Q))
        return tests[key]
    bits=[1<<i for i in range(len(orbits))]
    sizes=[sum(len(o) for b,o in zip(bits,orbits) if mask&b) for mask in range(1<<len(orbits))]
    candidates=sorted((m for m,v in enumerate(sizes) if 0<v<15), key=lambda m:-sizes[m])
    for mask in candidates:
        if any((mask & m)==mask for m,k,p in blockers): continue
        worst=row_endpoints(mask)
        part,k=map(int,np.unravel_index(worst.argmin(),worst.shape))
        table,required=exact_table(part,k)
        def impossible(m):
            return sum(max(row[i] for i,b in enumerate(bits) if b&m) for row in table)<required
        assert impossible(mask), (mask,part,k,float(worst[part,k]))
        # If this single angle still rules out a larger subset, it also rules out all of its subfamilies.
        for b in bits:
            if not mask&b and impossible(mask|b):mask|=b
        blockers.append((mask,k,part))
    assert all(any(mask & m == mask for m,k,p in blockers) for mask in candidates)
    out.mkdir(exist_ok=True,parents=True)
    (out/'fifteen.json').write_text(json.dumps(small,separators=(',',':'))+'\n')
    (out/'fifteen-rows.json').write_text(json.dumps(dict(cells=proposals),separators=(',',':'))+'\n')
    (out/'fifteen-unit.json').write_text(json.dumps(unit,separators=(',',':'))+'\\n')
    receipt=dict(native_table_contraction=str(scale),
                 contracted_diameter_squared=str(F(sqdiam,squashed["denominator"]**2)),
                 contracted_minimum_integer_gain=min(a-b for a,b in zip(shrunken[7::8],old_scores)),
                 exact_unit_sha256=sha(out/'fifteen-unit.json'),source_sha256=sha(bankdir/'new30.json'),old_points_sha256=sha(comp/'finite-point-leaf-input.json'),
                 old_cells_sha256=sha(comp/'finite-point-global-input.json'),
                 points=chosen, square_diameter_numerator=diam, denominator=small['denominator'],
                 cells=len(cells), rows=len(gains), minimum_native_gain=min(gains),maximum_native_gain=max(gains),scale=Q,
                 reflection_orbits=orbits,nonempty_symmetric_subsets_with_at_most_14_points=len(candidates),
                 blocking_supersets_and_counterexample_rows=blockers,minimum_in_this_subset_model=15,
                 exact_upper_projection_scale=scales,certification_precision=ctx.prec,
                 points_sha256=sha(out/'fifteen.json'),rows_sha256=sha(out/'fifteen-rows.json'))
    (out/'fifteen-proof.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(len(candidates),len(blockers),min(gains),max(gains),diam,flush=True)

if __name__=='__main__':
    run(*(Path(v) for v in sys.argv[1:]))
