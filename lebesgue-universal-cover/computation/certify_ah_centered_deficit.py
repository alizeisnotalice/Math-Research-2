"""Integrate the derivative of g(ell(x))-t(x) directly in the physical coordinates of the new cap."""
from fractions import Fraction as F
import contextlib,io
from certify_ah_chain_chord import geometry as seed_geometry
from certify_ah_margin_envelope import ANCHORED_CUTS,replace_line


def geometry(archive,extra,*,far,negative,eps):
    old,prefix,_=seed_geometry(archive,[(F('.427'),F('.00225'))],far=F('.116'),negative=F('.123'),eps=F('.0026'))
    pmax=max(p for p,_ in extra);assert F('.419')<pmax<F('.43') and (pmax*5000).denominator==(far*1000).denominator==(negative*5000).denominator==1
    prefix=replace_line(prefix,next(l for l in prefix.splitlines() if l.startswith('cuts=')),ANCHORED_CUTS[:-1]+''.join(f', (F("{p}"), F("{d}"))' for p,d in extra)+']')
    prefix=replace_line(prefix,next(l for l in prefix.splitlines() if l.startswith('wp, wf,')),f'wp, wf, wmax, eps = F("{pmax}"), F("{far}"),F("0.43"), F("{eps}")');prefix=replace_line(prefix,next(l for l in prefix.splitlines() if l.startswith('pbound = ')),f'pbound = F("{negative}")')
    mark='parameter_boxes=[]'; add='''graph_slopes=[]
for k in range(430):
    aa,bb=F(k,1000),F(k+1,1000)
    ell=[chart(add(z,mul(A(vv),D)))[0] for vv in (aa,bb)]
    ll=ell[0].union(ell[1]);pp=coord_hull(psi,ll)
    _,dl,dt=points(box(aa,bb));mm=dot(dt,D)/dot(dl,D)
    dp=psi_prime(add(l0,mul(ll,el)))
    graph_slopes.extend(((ll,mm),(pp,-mm/dp)))

def graph_slope(ell):
    assert ell>lmin and ell<lmax
    possible=[s for rr,s in graph_slopes if rr.overlaps(ell)];assert possible
    return min(s.lower() for s in possible).union(max(s.upper() for s in possible))
''';prefix=replace_line(prefix,mark,add+'\n'+mark)
    old='''    assert value > A("0.008"), (a, value)
    cell_deficit = max(zero,((own_height - old_height_lower(l)) / value.lower()).upper())''';new='''    slope=graph_slope(l)
    deficit0=graph_value(cl)-ct
    dw_def=slope*dot(gl,D)-dot(gt,D)
    dh_def=slope*dot(gl,normal)-dot(gt,normal)
    cell_deficit=max(zero,(deficit0+step*dw_def+ht*dh_def).upper())
    if value>0:
        cell_deficit=min(cell_deficit,max(zero,((own_height-old_height_lower(l))/value.lower()).upper()))''';prefix=replace_line(prefix,old,new)
    # l must contain the chart along both coordinate legs from (mid,0) to (w,h).
    ns={'__file__':str(archive/'certify_ah_two_line_cap.py'),'__name__':'__main__'}
    with contextlib.redirect_stdout(io.StringIO()):exec(compile(prefix,ns['__file__'],'exec'),ns)
    return ns,prefix,None
