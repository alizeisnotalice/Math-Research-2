"""Restrict each negative physical cap cell to the known paired arc's upper side."""
import contextlib
from fractions import Fraction as F
from functools import lru_cache
import io

from flint import ctx
from certify_ah_centered_deficit import geometry as centered_geometry
from certify_ah_margin_envelope import ANCHORED_CUTS, replace_line


@lru_cache(maxsize=2)
def template(archive, labelled):
    if labelled:
        from certify_ah_first_hit_chain import geometry as labelled_geometry
        ns, source, _ = labelled_geometry(
            archive, [(F('.427'), F('.0038'))], far=F('.16'), negative=F('.177'), eps=F('.0043'))
    else:
        ns, source, _ = centered_geometry(
            archive, [(F('.427'), F('.00225'))], far=F('.116'), negative=F('.123'), eps=F('.0026'))
    return ns['__file__'], source


FLOOR = '''physical_boundaries=[]
physical_lowers=[]
for k in range(int(5000*wf)+1):
    t=A(F(k,5000))
    physical_boundaries.append(dot(direct,sub(q(t),z))/length)
assert all(x<y for x,y in zip(physical_boundaries,physical_boundaries[1:]))
for k in range(int(5000*wf)):
    a,b=F(k,5000),F(k+1,5000)
    _,derivatives=point_jets(box(a,b))
    dq=reflection(derivatives[2])
    assert dot(direct,dq)>0
    at=dot(normal,sub(q(A((a+b)/2)),z))
    physical_lowers.append(max(zero,(at+box((a-b)/2,(b-a)/2)*dot(normal,dq)).lower()))
physical_lefts=[x.lower() for x in physical_boundaries]
physical_rights=[x.upper() for x in physical_boundaries]

def physical_floor(a,b):
    if a>=0:
        return zero
    if b>=0:
        return zero
    start=max(0,bisect_left(physical_rights,A(-b).lower())-1)
    finish=min(len(physical_lowers),bisect_right(physical_lefts,A(-a).upper())+1)
    assert start<finish
    return min(physical_lowers[start:finish])

'''


def geometry(archive, extra, *, far, negative, eps, labelled=False, two_hit=False):
    ctx.prec=128
    path, source = template(archive,labelled)
    right = max(p for p,_ in extra)
    assert F('.419')<right<F('.43') and right*5000==int(right*5000)
    assert far*1000==int(far*1000) and negative*5000==int(negative*5000)
    cuts=ANCHORED_CUTS[:-1]+''.join(f', (F("{p}"), F("{d}"))' for p,d in extra)+']'
    source=replace_line(source,next(l for l in source.splitlines() if l.startswith('cuts=')),cuts)
    source=replace_line(source,next(l for l in source.splitlines() if l.startswith('wp, wf,')),
                        f'wp, wf, wmax, eps = F("{right}"), F("{far}"), F(".43"), F("{eps}")')
    source=replace_line(source,next(l for l in source.splitlines() if l.startswith('pbound = ')),
                        f'pbound = F("{negative}")')
    if two_hit:
        tube='tuple((eps*F(i,2),eps*F(i+1,2)) for i in range(-2,2))'
        assert source.count(tube)==2
        source=source.replace(tube,'tuple((eps*F(i,6),eps*F(i+1,6)) for i in range(-6,6))')
    source=replace_line(source,'parameter_boxes=[]',FLOOR+'parameter_boxes=[]')
    before='''    own_height = height_bound(a)
    x=add(add(z,mul(w,D)),mul(box(0,own_height.upper()),normal))
    cl,ct=chart(add(z,mul(A((a+b)/2),D)))'''
    after='''    if A(-b)>physical_rights[-1]:
        continue
    own_height=height_bound(a)
    floor=physical_floor(a,b)
    if own_height<floor:
        continue
    x=add(add(z,mul(w,D)),mul(floor.union(own_height.upper()),normal))
    cl,ct=chart(add(add(z,mul(A((a+b)/2),D)),mul(floor,normal)))'''
    source=replace_line(source,before,after)
    source=replace_line(source,'ht=box(0,own_height.upper())','ht=box(0,(own_height-floor).upper())')
    if two_hit:
        assert labelled
        source=replace_line(source,'parameter_boxes=[]','physical_deficits=[]\nparameter_boxes=[]')
        marker='    parameter_boxes.append((l,pl,used,a))'
        source=replace_line(source,marker,'    physical_deficits.append((l,cell_deficit))\n'+marker)
    ns={'__file__':path,'__name__':'__main__'}
    with contextlib.redirect_stdout(io.StringIO()):
        exec(compile(source,path,'exec'),ns)
    if labelled:
        assert all(ns['deficit_by_bin'][j]==max(row[j] for row in ns['labelled_deficits'])
                   for j in range(len(ns['ref_boxes'])))
    if two_hit:
        ns['two_hit_budgets']={}
        A=ns['A']
        for j,E in enumerate(ns['deficit_by_bin']):
            if E.is_zero():
                continue
            a,b=ns['knots'][j:j+2]
            l1,l2=map(lambda w: ns['chart'](ns['add'](ns['z'],ns['mul'](A(w),ns['D'])))[0],(a,b))
            ell=l1.union(l2)
            partner=ns['coord_hull'](ns['psi'],ell)
            budgets=[]
            for q,center in enumerate((partner,ell)):
                first=ns['labelled_deficits'][q][j]
                locus=center+(-A(F(2,3))*first).union(A(0))
                second=max([A(0)]+[d for x,d in ns['physical_deficits'] if x.overlaps(locus)])
                budgets.append((first,second))
            ns['two_hit_budgets'][j]=budgets
    return ns,source,None
