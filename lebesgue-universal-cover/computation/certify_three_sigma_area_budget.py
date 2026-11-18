"""Rebuild the same AH and eleven EH cuts at three exact Gibbs angles."""
from argparse import ArgumentParser
from fractions import Fraction as F
from hashlib import sha256
import json
import os
from itertools import combinations
from pathlib import Path
import subprocess
import sys

from flint import arb,ctx

ANGLES=('1.54','1.5497','1.56')
BASELINE='1.5497'
PIVOT=F('.427')
TURN=F('.006')
U0=F('0.844092046521218581802948589770767704354')


def digest(p):
    return sha256(p.read_bytes()).hexdigest()


def exact(x):
    f=F(x)
    return arb(f.numerator)/f.denominator


def rational(x):
    a,b=map(int,x.man_exp())
    return F(a)*F(2)**b


def enclosure(x):
    return rational(x.lower()),rational(x.upper())


def pack(x,places=22):
    if isinstance(x,arb):
        x=enclosure(x)
    lo,hi=x
    k=10**places
    limits=(int((lo*k).__floor__()),int((hi*k).__ceil__()))
    def decimal(a):
        sign='-' if a<0 else ''
        z=abs(a)
        return f'{sign}{z//k}.{z%k:0{places}d}'
    return {'decimal':[decimal(z) for z in limits],'exact':[str(lo),str(hi)]}


def subtract(a,b):
    return a[0]-b[1],a[1]-b[0]


def specialize(source,token,value):
    assert source.count(token)==1,token
    return source.replace(token,value)


def fixed_ah(archive,sigma,pivot=PIVOT,turn=TURN,upper=F(".203")):
    source=archive/'certify_ah_combined_area.py'
    code=source.read_text().split('baseline = here / ',1)[0]
    for a,b in (
        ('sigma = "1.5497"',f'sigma = "{sigma}"'),
        ('pivot_parameter = F("0.38")',f'pivot_parameter = F("{pivot}")'),
        ('angle = F("0.00033")',f'angle = F("{turn}")'),
        ('start, finish = F("0.037"), F("0.038")',f'start, finish = F(0), F("{upper}")'),
    ):
        code=specialize(code,a,b)
    ns={'__file__':str(source),'__name__':'__main__'}
    exec(compile(code,str(source),'exec'),ns)
    z,D,W,ld=[ns[k] for k in ('z','D','W','line_direction')]
    assert ns['cross'](ld,ns['sub'](z,W)).hi<0
    assert ns['cross'](ld,ns['sub'](ns['add'](z,D),W)).lo>0
    assert ns['cut'](F(0)).hi<0 and ns['cut'](F(1)).lo>0
    thin=ns['cap_area']
    return (thin.lo,thin.hi),(ns['start'],ns['finish']),ns['slope'].lo


def signed_arc(a,b,c,g):
    u,v=g.minus(a,c),g.minus(b,c)
    dot=g.dot(u,v)
    assert dot>0
    return g.cross(c,g.minus(b,a))+(g.cross(u,v)/dot).atan()


def setup(sigma,g):
    g.PI=arb.pi()
    g.SIGMA=exact(sigma)*g.PI/180
    n,m,C3,B3,_,U,_,T,J=g.setup()
    meet,rot,sl,dot=g.line_pair,g.rot,g.scale,g.dot
    sub,plus,cl,cc=g.minus,g.plus,g.circle_line,g.circle_pair
    def project(a,w):
        return sub(a,sl(dot(a,w)-g.HALF,w))
    H=project(C3,n[3])
    Ip=project(B3,m[3])
    E3,E2,C2,F3=(meet(n[i],m[j]) for i,j in ((4,3),(3,3),(5,5),(3,2)))
    A1=meet(n[1],n[2])
    K=project(F3,m[5])
    G=cl(F3,n[5],1)
    Q,P=project(E3,n[1]),project(G,n[2])
    X=cc(G,E3,-1)
    c=(-E3[0],E3[1])
    r3=arb(3).sqrt()
    z=cc((r3/2-g.HALF,(r3/2-g.HALF-1)/r3),(-g.HALF/2,-r3/4),-1)
    zr=sub(z,c)
    zz=dot(zr,zr)
    Ta=plus(c,sl(1/zz,sub(zr,sl((zz-1).sqrt(),rot(zr)))))
    return locals()


def ah_clearance(V,upper,T,g,pivot,turn):
    z,Ta=V['z'],V['Ta']
    delta=g.minus(Ta,z)
    b=g.scale(-1/g.norm(delta),delta)
    c,s=g.unit(exact(turn))
    cut=g.plus(g.scale(c,b),g.scale(s,g.rot(b)))
    anchor=g.plus(z,g.scale(exact(pivot),delta))
    assert g.cross(cut,g.minus(z,anchor))<0 and g.cross(cut,g.minus(Ta,anchor))>0
    margins=[g.cross(cut,g.minus((x,y),anchor)) for x in (upper[0],T[0]) for y in (T[1],upper[1])]
    assert all(v>0 for v in margins)
    return min(v.lower() for v in margins)


def geometries(sigma,old,parameters,g):
    V=setup(sigma,g)
    C3,B3,U,T,J=[V[k] for k in ('C3','B3','U','T','J')]
    plus,minus,scale,cross,dot,rot=g.plus,g.minus,g.scale,g.cross,g.dot,g.rot
    def arc(a,b,c):
        return signed_arc(a,b,c,g)
    def moved(d):
        rad=minus(U,C3)
        c,s=g.unit(exact(d))
        return plus(C3,minus(scale(c,rad),scale(s,rot(rad))))
    def halfplane(a,b):
        e=minus(b,a)
        u=scale(-1/g.norm(e),rot(e))
        return u,dot(u,a)
    def intersect(row1,row2):
        return g.line_pair(row1[0],row2[0],row1[1],row2[1])
    def ccw_window(a,b,c):
        aa,bb=minus(a,c),minus(b,c)
        assert dot(aa,bb)>0 and cross(aa,bb)>0
    H,Ip=V['H'],V['Ip']
    d0,d1,h=map(F,old)
    u,u1=moved(d0),moved(d1)
    hinge=plus(u1,scale(exact(h),minus(T,u1)))
    assert cross(minus(hinge,u),minus(T,hinge))>0
    r126=[halfplane(u,hinge),halfplane(hinge,T)]
    from certify_eh_dual_hinge import RELAX
    relax=RELAX
    X=intersect(*[(n,c-exact(d)) for (n,c),d in zip(r126,relax,strict=True)])
    rows=[halfplane(u,X),halfplane(X,T)]
    upper_corners=[]
    last=d0
    prev_p=F(0)
    for d,p in map(lambda r:map(F,r),parameters):
        assert d>last and prev_p<p<1
        earlier=moved(d)
        Y=plus(X,scale(exact(p),minus(u,X)))
        rows.append(halfplane(earlier,Y))
        upper_corners.append(earlier)
        last,prev_p=d,p
    assert len(rows)==11
    upper=upper_corners[-1]
    ccw_window(H,upper,C3)
    ccw_window(upper,U,C3)
    ccw_window(U,J,C3)
    ccw_window(J,T,B3)
    ccw_window(T,Ip,B3)
    assert all((q[0]<0 and q[1]<0) for q in [minus(upper,C3),minus(U,C3)])
    assert upper[0]<T[0] and upper[1]>T[1]
    def slacks(z,exclude):
        return [b-dot(n,z) for i,(n,b) in enumerate(rows) if i not in exclude]
    vertices=[(upper,{10}),(T,{1})]
    for (i,r),(j,s) in combinations(enumerate(rows),2):
        z=intersect(r,s)
        if z[1]<T[1] or z[1]>upper[1]:continue
        assert T[1]<z[1]<upper[1]
        residual=slacks(z,{i,j})
        assert all(a>0 or a<0 for a in residual)
        if all(a>0 for a in residual):
            vertices.append((z,{i,j}))
    vertices.sort(key=lambda r:r[0][1].mid(),reverse=True)
    end_slacks=slacks(upper,{10})+slacks(T,{1})
    assert all(s>0 for s in end_slacks)
    points=[z for z,_ in vertices]
    edges=[minus(b,a) for a,b in zip(points,points[1:])]
    facets=[]
    for (a,I),(b,J0) in zip(vertices,vertices[1:]):
        common=I&J0
        assert len(common)==1 and b[1]<a[1]
        facets.extend(common)
    assert len(facets)==len(set(facets))
    assert all(a[0]>0 and a[1]<0 for a in edges)
    assert all(cross(a,b)>0 for a,b in zip(edges,edges[1:]))
    assert cross(rot(minus(upper,C3)),edges[0])>0
    assert cross(edges[-1],rot(minus(T,B3)))>0

    # Every intermediate vertex is in the triangle of three points on the
    # original convex body. This certifies the inserted chain belongs to it.
    triangle=[upper,U,T]
    for a in points[1:-1]:
        tests=[cross(minus(r,s),minus(a,s)) for s,r in zip(triangle,triangle[1:]+triangle[:1])]
        assert all(t>0 for t in tests),[str(t) for t in tests]

    clearance=ah_clearance(V,upper,T,g,PIVOT,TURN)
    old_eh=(cross(T,U)+arc(U,J,C3)+arc(J,T,B3))/2
    new_eh=(arc(upper,J,C3)+arc(J,T,B3)-
            sum(cross(a,b) for a,b in zip(points,points[1:])))/2
    r126_eh=(arc(u,J,C3)+arc(J,T,B3)-cross(u,hinge)-cross(hinge,T))/2
    assert new_eh>r126_eh>old_eh>0
    def lines(*names):
        pts=[V[k] for k in names]
        return sum(cross(a,b) for a,b in zip(pts,pts[1:]))
    cs=(lines('G','C2','K')+arc(V['K'],V['G'],V['F3']))/2
    es=(lines('H','E2','Ip')+arc(V['Ip'],J,B3)+arc(J,H,C3))/2
    a_s=(lines('P','A1','Q')+arc(V['Q'],V['X'],V['E3'])+arc(V['X'],V['P'],V['G']))/2
    cos=g.SIGMA.cos()
    p_area=arb(3).sqrt()/2-2*arb(3).sqrt()*(cos/arb(3).sqrt()-g.HALF)**2/(4*cos**2-3)
    assert all(v>0 or v<0 for v in (cs,es,a_s))
    cs,es,a_s=map(abs,(cs,es,a_s))
    return dict(p_area=p_area,cs=cs,es=es,a_s=a_s,old_eh=old_eh,new_eh=new_eh,
                r126_eh=r126_eh,facets=facets,clearance=clearance,
                points=points)


def run(a):
    sys.path.insert(0,str(a.comp.resolve()))
    import audit_optimized_eh_chord as g
    ctx.prec=192
    archive=a.comp
    read=lambda n:json.loads((archive/n).read_bytes())
    old=read('eh-two-line-cover-certificate.json')
    parameters=read('eh-upper-chain-parameters.json')
    r199=read('flexible-006-eleven-eh.json')
    assert old['parameters']==['0.000555','0.0005522','0.0075']
    assert len(parameters)==9
    r199_ah=read('ah-flexible-chord-certificate.json')
    assert [[F(x) for x in row] for row in r199_ah['lines']]==[[PIVOT,TURN]]
    assert r199['AH_sha256']==digest(archive/'ah-flexible-chord-certificate.json')
    assert r199['parameter_sha256']==digest(archive/'eh-upper-chain-parameters.json')
    assert r199['R126_sha256']==digest(archive/'eh-two-line-cover-certificate.json')
    base_area_ref=tuple(map(F,r199['area_exact_interval']))
    summaries=[]
    extra_rows=[]
    for sigma in ANGLES:
        env={**os.environ,'LEBESGUE_SIGMA_DEGREES':f'{sigma},{sigma}'}
        sp=archive/'certify_splice_connectors.py'
        shape=json.loads(subprocess.check_output([sys.executable,'-B',str(sp)],env=env))
        assert shape['passed'] and shape['sigma_degrees_interval']==[str(F(sigma))]*2
        basepy=archive/'certify_gibbs_area_refined.py'
        base=json.loads(subprocess.check_output([sys.executable,'-B',str(basepy),sigma,'8','32']))
        assert base['passed'] and F(base['sigma_degrees'])==F(sigma)
        gibbs_area=tuple(map(F,base['area_interval']))
        geo=geometries(sigma,old['parameters'],parameters,g)
        ah,w,slope=fixed_ah(archive,sigma)
        remainder=geo['p_area']-geo['cs']-geo['es']-geo['a_s']-geo['old_eh']
        old_ah=subtract(enclosure(remainder),gibbs_area)
        assert old_ah[0]>0
        extra_eh=geo['new_eh']-geo['old_eh']
        gross_eh=geo['new_eh']
        area=subtract(subtract(gibbs_area,ah),enclosure(extra_eh))
        if sigma==BASELINE:
            assert area[0]<base_area_ref[0]<=base_area_ref[1]<area[1]
            assert geo['facets']==r199['active_facet_indices']
            known_AH=tuple(map(F,r199_ah['shape']['cap_area_exact_interval']))
            known_EH=tuple(map(F,r199['EH_gain_over_R126_exact_interval']))
            known_126=tuple(map(F,old['cut_gain_area_interval_exact']))
            for fresh,formerly in (
                (ah,known_AH),
                (enclosure(geo['new_eh']-geo['r126_eh']),known_EH),
                (enclosure(geo['r126_eh']-geo['old_eh']),known_126)
            ):
                assert max(fresh[0],formerly[0])<=min(fresh[1],formerly[1])

        summaries.append(dict(sigma=sigma,base=pack(geo['p_area']),
            sprague={k:pack(geo[k]) for k in ('cs','es','a_s')},
            old_ah=pack(old_ah),old_eh=pack(geo['old_eh']),
            r34_gibbs=pack(gibbs_area),extra_ah=pack(ah),
            extra_eh=pack(extra_eh),absolute_eh=pack(gross_eh),
            extra_eh_over_r126=pack(geo['new_eh']-geo['r126_eh']),
            eh_active=geo['facets'],eh_boundary_vertices=[[pack(q) for q in p] for p in geo['points']],
            eh_clearance_to_ah=pack(geo['clearance']),ah_crossing=pack(w),
            ah_crossing_slope_lower=str(slope),combined=pack(area),
            saving_from_r199=pack(subtract(base_area_ref,area)),
            more_cut_needed_for_A=pack(subtract(area,(U0-F('1e-6'),)*2)),
            more_cut_needed_for_B=pack(subtract(area,(U0-F('1e-5'),)*2))))
        if sigma=='1.54' and a.bracket_a:
            pivot=F('.4285')
            upper,T=geo['points'][0],geo['points'][-1]
            V=setup(sigma,g)
            for turn in map(F,('.006','.00890','.00893')):
                cap,root,slope=fixed_ah(archive,sigma,pivot,turn,F('.30'))
                alt_area=subtract(subtract(gibbs_area,cap),enclosure(extra_eh))
                margin=ah_clearance(V,upper,T,g,pivot,turn)
                extra_rows.append((turn,alt_area,dict(
                    pivot=str(pivot),turn=str(turn),AH_cap=pack(cap),
                    crossing=pack(root),AH_EH_clearance=pack(margin),area=pack(alt_area),
                    saving_from_r199=pack(subtract(base_area_ref,alt_area)),
                    saving_from_A_target=pack(subtract((U0-F('1e-6'),)*2,alt_area)))))
    if a.bracket_a:
        assert len(extra_rows)==3
        reference,lo,hi=[z[1] for z in extra_rows]
        assert reference[0]>base_area_ref[1]
        assert lo[0]>U0-F('1e-6')>hi[1]
    return dict(cases=summaries,area_only_A_bracket=[z[2] for z in extra_rows],r199_area=pack(base_area_ref),pivot=str(PIVOT),turn=str(TURN),
                eh_original_parameters=old['parameters'],eh_nine_additional=parameters,
                precision=ctx.prec,source_sha256=digest(Path(__file__)),
                sources={n:digest(archive/n) for n in (
                    'certify_splice_connectors.py','certify_ah_combined_area.py','certify_gibbs_area_refined.py',
                    'audit_optimized_eh_chord.py','certify_eh_dual_hinge.py','eh-two-line-cover-certificate.json',
                    'eh-upper-chain-parameters.json','ah-flexible-chord-certificate.json','flexible-006-eleven-eh.json')})


if __name__=='__main__':
    assert __debug__
    parser=ArgumentParser(description=__doc__)
    parser.add_argument('comp',type=Path)
    parser.add_argument('--bracket-a',action='store_true')
    print(json.dumps(run(parser.parse_args()),separators=(',',':')))
