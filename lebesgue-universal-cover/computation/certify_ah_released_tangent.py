"""A convex two-tangent replacement for the AH paired arc. Prove true area and full-body placements."""
import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
import sys

from flint import arb, arb_series, ctx, __version__ as flint_version


def load_geometry(archive):
    source = archive / 'certify_ah_two_line_cap.py'
    code = source.read_text()
    ns = {'__file__': str(source)}
    # Reuse the archived five vertices, exact cone ordering and support formulas.
    prefix = code.split('for pivot,c in zip(pivot_set,cut_dirs):')[0]
    exec(compile(prefix, str(source), 'exec'), ns)
    ns['sector'] = ns['audit']['sector_formula']
    ns['orders'] = ns['audit']['order']
    exec(compile(code[code.index('def rotate('):code.index('def F_at(')], str(source), 'exec'), ns)
    return ns


def triangle_support(g,A,u):
    r=1/A(3).sqrt()
    vv=[(A(0),r),(A('.5'),-r/2),(-A('.5'),-r/2)]
    definite=[g.dot(x,u) for x in vv]
    possible=definite.copy()
    for i,x in enumerate(vv):
        d,e=[g.minus(y,x) for j,y in enumerate(vv) if j!=i]
        if g.cross(d,e)<0:
            d,e=e,d
        assert g.cross(d,e)>0
        left,right=g.cross(d,u),g.cross(u,e)
        radial=1+g.dot(x,u)
        if left>=0 and right>=0:
            definite.append(radial)
        if not (left<0 or right<0):
            possible.append(radial)
    return max(y.lower() for y in definite).union(max(y.upper() for y in possible))


def certify(archive, p, delta, i, j, bounds, tmax, inner_offset=F(0)):
    import audit_optimized_eh_chord as g
    from audit_eh_continuous_placement import polygon
    from certify_eh_upper_chain import front
    from certify_eh_joint_parameters import specialize
    ns = load_geometry(archive)
    A = ns['A']
    specialize(F('1.5497'))
    ctx.prec = 160
    ctx.cap = 10
    def read(name):
        return json.loads((archive/name).read_bytes())
    ndiv = 128
    a,b = F(i,ndiv),F(j,ndiv)
    oldp, oldangle = F('.427'), F('.006')
    assert 0 < oldp < p <= 1 and 0 < oldangle < delta < 1
    t = A(tmax)
    def q(w):
        return ns['reflection'](ns['points'](w)[0][2])
    def tangent(w):
        pair = q(arb_series([A(w),1]))
        point = tuple(x[0] for x in pair)
        der = tuple(x[1] for x in pair)
        normal = g.scale(-1/g.norm(der),g.rot(der))
        return point, normal, g.dot(point,normal)
    left, ul, hl = tangent(a)
    right, ur, hr = tangent(b)
    tangent_rows=[('left_tangent',ul,hl)]
    if inner_offset:
        for k in range(i+4,j-3,4):
            _,u,h=tangent(F(k,ndiv))
            tangent_rows.append((f'offset_tangent_{k}',u,h+A(inner_offset)))
    tangent_rows.append(('right_tangent',ur,hr))
    corners=[g.line_pair(u,v,h,beta) for (_,u,h),(_,v,beta) in zip(tangent_rows,tangent_rows[1:])]
    corner=corners[0]
    allpts=[left,*corners,right]
    hexagon, tilted, *_=g.setup()
    Pslack=min((A('.5')-g.dot(normal,x)).lower() for x in allpts for normal in (*hexagon,tilted[3],tilted[5]))
    assert Pslack>0
    for x,y,z in zip(allpts,allpts[1:],allpts[2:]):
        assert g.cross(g.minus(y,x),g.minus(z,y))>0
    assert all(g.cross(u,v)>0 for (_,u,_),(_,v,_) in zip(tangent_rows,tangent_rows[1:]))
    assert all(g.dot(g.minus(y,x),g.rot(u))>0 for x,y,(_,u,_) in zip(allpts,allpts[1:],tangent_rows))

    assert g.cross(ul,ur)>0
    assert g.cross(g.minus(corner,left),g.minus(right,corner))>0
    direction = ns['direct']
    angle = A(delta)
    cutnormal = g.neg(g.rot(g.plus(g.scale(angle.cos(),direction),g.scale(angle.sin(),g.rot(direction)))))
    W = g.plus(ns['z'],g.scale(A(p),ns['D']))
    cutheight = g.dot(cutnormal,W)
    def cut_slack(x):
        return cutheight-g.dot(cutnormal,x)
    # The outward patch is inside this triangle. The two clipping operations are disjoint.
    patch_slack = min(cut_slack(x).lower() for x in allpts)
    assert patch_slack > 0
    lo,hi = F(0),a
    assert cut_slack(q(A(lo)))<0<cut_slack(q(A(hi)))
    for _ in range(105):
        mid = (lo+hi)/2
        score = cut_slack(q(A(mid)))
        if score < 0:
            lo=mid
        else:
            assert score>0
            hi=mid

    def density(pos,base):
        x,y=q(arb_series([pos,1]))
        return (x-base[0])*y.derivative()-(y-base[1])*x.derivative()
    def green(l,r,base):
        result = A(0)
        # 8 even-order quadrature with the 8th Taylor coefficient enclosed on its whole cell.
        for k in range(128):
            u=l+(r-l)*k/128
            v=l+(r-l)*(k+1)/128
            half=A((v-u)/2)
            span=A(u).union(A(v))
            center=density(A((u+v)/2),base)
            rem=density(span,base)[8]
            assert rem.is_finite()
            result += sum(center[d]*half**(d+1)/(d+1) for d in range(0,8,2))
            result += rem*half**9/9
        return result
    cap = green(F(0),lo,W)+A(0).union(A((hi-lo)/2))*density(A(lo).union(A(hi)),W)[0]
    added=sum(g.cross(g.minus(x,left),g.minus(y,left))/2 for x,y in zip(allpts,allpts[1:]))-green(a,b,left)
    assert cap>0 and added>0
    rows=polygon('0.0005522',ndiv)[0]
    rows=[(name,u,h) for name,u,h in rows if not (name.startswith('chord_Z_') and i<=int(name.split('_')[-1])<j)]
    rows.extend(tangent_rows)
    rows.append(('new_AH',cutnormal,cutheight))
    old=read('eh-two-line-cover-certificate.json')
    eh=front(old,read('eh-upper-chain-parameters.json'))
    rows.extend((f'EH_{k}',u,h) for k,(u,h) in enumerate(zip(eh['normals'],eh['offsets'])))
    # Verify that the exact unit-distance formulas have the same ten normal cones
    # and the same three hexagon contacts over the entire proposed w interval.
    for cell in range(512):
        w0=bounds[0]+(bounds[1]-bounds[0])*cell/512
        w1=bounds[0]+(bounds[1]-bounds[0])*(cell+1)/512
        v=ns['centered_points'](w0,w1)
        edges=[g.minus(v[k],v[l]) for k,l in ns['orders']]
        for k,e in enumerate(edges):
            f=edges[(k+1)%10]
            assert (e[1]>0 if k<5 else e[1]<0)
            assert g.cross(e,f)>0 and g.dot(e,f)>A('.75')
        for k in range(5):
            for l in range(k+1,5):
                if (k,l) not in ns['orders'] and (l,k) not in ns['orders']:
                    sep=g.minus(v[k],v[l])
                    assert g.dot(sep,sep)<1
        for k,u in enumerate(hexagon):
            cr=[g.cross(e,u) for e in edges]
            possible={ns['sector'][j] for j in range(10) if not cr[j]<0 and not cr[(j+1)%10]>0}
            assert possible=={((1,3,0)[k%3],int(k<3))}
    # The six unrotated hexagon strips are exact equalities of the normalized constant-width bodies.
    exact_hex={'B1-Q','P-F1','F1-H','E3-D1','D1-G','C3-B1'}
    hexrows=[(name,u,h) for name,u,h in rows if name in exact_hex]
    assert len(hexrows)==6
    for name,u,h in hexrows:
        assert (h-A('.5')).contains(0)
        assert sum(all((x-y).contains(0) for x,y in zip(u,vv)) for vv in g.setup()[0])==1
    active=[x for x in rows if x[0] not in exact_hex and not (inner_offset and x[0].startswith(('chord_Tp_','chord_G_','chord_Te_','chord_H_','chord_Q_')))]
    assert {'I-E3','K-C3'}<=set(r[0] for r in active)
    # Triangle supports are evaluated on their three circle fans. Minkowski linearity
    # makes the two endpoint support bounds valid throughout the closed weight interval.
    R=1/A(3).sqrt()
    def slacks(w0,w1,rr):
        mid=A((w0+w1)/2)
        v=ns['centered_points'](w0,w1)
        vv=ns['points'](mid)[0]
        dv=ns['point_jets'](A(w0).union(A(w1)))[1]
        rad=A((w1-w0)/2)
        edges=[g.minus(v[k],v[l]) for k,l in ns['orders']]
        ans=[]
        for name,u,h in rr:
            pre=ns['reflection'](u)
            cr=[g.cross(e,pre) for e in edges]
            feasible=[(g.dot(vv[k],pre)+c+rad*abs(g.dot(dv[k],pre))).upper()
                      for z,(k,c) in enumerate(ns['sector']) if not cr[z]<0 and not cr[(z+1)%10]>0]
            assert feasible
            val=max(feasible)
            tr=triangle_support(g,A,u)
            hi=val+t*(tr-val).union(A(0))
            ans.append((h-hi).lower())
        return ans
    # First discharge all easy rows on each coarse interval, then refine only unresolved rows.
    cells=0
    worst=A(1)
    def sweep(l,r,rr):
        nonlocal cells,worst
        cells+=1
        candidates=[]
        for row,s in zip(rr,slacks(l,r,rr)):
            if s>0:
                worst=min(worst,s)
            else:
                candidates.append(row)
        if candidates:
            assert r-l>F(1,10**9),(l,r,[s[0] for s in candidates])
            m=(l+r)/2
            sweep(l,m,candidates)
            sweep(m,r,candidates)
    sweep(*bounds,active)
    vals=slacks(F('.429'),F('.429'),active)
    explicit=dict((row[0],str(s)) for row,s in zip(active,vals) if row[0] in ('left_tangent','right_tangent','new_AH'))
    assert min(vals)>0
    oldnormal = tangent(F('.429'))[1]
    V=(A(0),R)
    oldexcess=A('1e-8')*g.dot(oldnormal,g.minus(V,q(A('.429'))))
    assert oldexcess>A('9e-11')
    fits={}
    for filename,n in [('widened-projected-obstruction-points.json',60),('widened-projected-obstruction-18-points.json',18)]:
        record=read(filename)
        pts=[tuple(A(F(z,record['denominator'])) for z in xy) for xy in record['points']]
        assert len(pts)==n
        minimum=min((h-g.dot(u,x)).lower() for name,u,h in rows for x in pts)
        fits[filename]=str(minimum)
    assert all(arb(s)>0 for s in fits.values())
    base=read('flexible-006-eleven-eh.json')
    oldah=read('ah-flexible-chord-certificate.json')
    assert base['AH_sha256']==sha256((archive/'ah-flexible-chord-certificate.json').read_bytes()).hexdigest()
    assert [list(map(F,x)) for x in oldah['lines']]==[[oldp,oldangle]]
    enc=lambda v:A(F(v[0])).union(A(F(v[1])))
    oldarea=enc(base['area_exact_interval'])
    oldcap=enc(oldah['shape']['cap_area_exact_interval'])
    deleted=cap-oldcap
    area=oldarea-deleted+added
    best=enc(read('fan-0062-eleven-eh.json')['area_exact_interval'])
    def fmt(x):
        # Exact outward decimal endpoints, including the uncertainty in all stored interval data.
        scale=10**35
        return [str(F(int((endpoint*scale).floor().unique_fmpz() if lower else (endpoint*scale).ceil().unique_fmpz()),scale))
                for lower,endpoint in ((True,x.lower()),(False,x.upper()))]
    return dict(p=str(p),angle_radians=str(delta),patch=[str(a),str(b)],interior_tangent_offset=str(inner_offset),patch_tangent_rows=len(tangent_rows),fold_family_w=list(map(str,bounds)),
                perturbation_weight_closed=['0',str(tmax)],subdivisions=ndiv,precision=160,
                patch_deep_line_slack=str(patch_slack),paired_crossing=list(map(str,(lo,hi))),
                added_area=str(added),deleted_area=str(deleted),overlap_correction_exact='0',candidate_true_area=str(area),
                true_area_rational=fmt(area),gain_against_R199=str(oldarea-area),gain_against_0062=str(best-area),
                explicit_full_body_slacks=explicit,explicit_old_tangent_excess=str(oldexcess),
                uniform_full_body_lower_slack=str(worst),verified_cells=cells,normal_fan_cells=512,hex_equalities=sorted(exact_hex),
                all_other_rows=len(active),new_boundary_P_clearance=str(Pslack),historical_identity_slacks=fits,
                true_Sprague_contact_groups=['Q','Tp','H','Te','G'] if inner_offset else [],universal_placement_proved=False)


if __name__=='__main__':
    cli=argparse.ArgumentParser()
    cli.add_argument('archive',type=Path)
    cli.add_argument('output',type=Path)
    cli.add_argument('--repair',action='store_true')
    cli.add_argument('--long',action='store_true')
    opt=cli.parse_args()
    archive=opt.archive.resolve()
    sys.path.insert(0,str(archive))
    if opt.long:
        result=certify(archive,F('.950'),F('.0110'),53,128,(F('.4285'),F('.950')),F('.00005'),F('.000001'))
    elif opt.repair:
        result=certify(archive,F('.430'),F('.0175'),53,58,(F('.4285'),F('.4301')),F('.00005'))
    else:
        result=certify(archive,F('.430'),F('.018'),54,56,(F('.42868'),F('.430')),F('1e-8'))
    names=['certify_ah_two_line_cap.py','ah_root_geometry_at_angle.py','audit_eh_continuous_placement.py',
           'audit_optimized_eh_chord.py','certify_eh_upper_chain.py','eh-upper-chain-parameters.json',
           'eh-two-line-cover-certificate.json','flexible-006-eleven-eh.json','fan-0062-eleven-eh.json',
           'ah-flexible-chord-certificate.json','widened-projected-obstruction-points.json',
           'widened-projected-obstruction-18-points.json']
    result['source_sha256']={name:sha256((archive/name).read_bytes()).hexdigest() for name in names}
    result['verifier_sha256']=sha256(Path(__file__).read_bytes()).hexdigest()
    result['python']=sys.version
    result['flint']=flint_version
    opt.output.write_text(json.dumps(result,indent=2)+'\n')
    print(result['candidate_true_area'],result['gain_against_R199'],result['verified_cells'])
