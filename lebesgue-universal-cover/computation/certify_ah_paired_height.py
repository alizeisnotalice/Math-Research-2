"""Subtract a rigorous lower height of the paired graph before assigning chart deficits."""
from fractions import Fraction as F
from pathlib import Path
import contextlib,io,inspect
from certify_ah_margin_envelope import ANCHORED_CUTS, replace_line
from certify_ah_wide_contact_chart import geometry as wide_geometry


def geometry(archive,extra,*,far=F(".10"),negative=F(".105"),eps=F(".0026")):
    # Build the verified wide source. Its graph, geometry and root references are unchanged.
    seed=((F('.427'),F('.00135')),)
    old,full=wide_geometry(archive,seed)
    top=max(p for p,_ in extra); assert F('.419')<top<F('.430') and (top*5000).denominator==1
    assert (far*1000).denominator==1 and (negative*5000).denominator==1
    orig=next(line for line in full.splitlines() if line.startswith('cuts='))
    cuts=ANCHORED_CUTS[:-1]+''.join(f', (F("{p}"), F("{d}"))' for p,d in extra)+']'
    full=replace_line(full,orig,cuts)
    initial=next(line for line in full.splitlines() if line.startswith('wp, wf,'))
    full=replace_line(full,initial,f'wp, wf, wmax, eps = F("{top}"), F("{far}"), F("0.43"), F("{eps}")')
    initial=next(line for line in full.splitlines() if line.startswith('pbound = ')); full=replace_line(full,initial,f'pbound = F("{negative}")')
    prefix=full.split('angular_cells=32')[0]
    marker='parameter_boxes=[]'
    height_table='''from bisect import bisect_left,bisect_right
paired_parameters=[]
paired_height=[]
for k in range(2151):
    e,t=chart(add(z,mul(A(F(k,5000)),D)))
    paired_parameters.append(psi(e))
assert all(a<b for a,b in zip(paired_parameters,paired_parameters[1:]))
for k in range(2150):
    a,b=F(k,5000),F(k+1,5000)
    midw=F(2*k+1,10000)
    middle=dot(normal,sub(q(midw),z))
    _,jets=point_jets(box(a,b))
    slope=dot(normal,reflection(jets[2]))
    enclosure=middle+box((a-b)/2,(b-a)/2)*slope
    paired_height.append(max(zero,enclosure.lower()))
low_ell=[x.lower() for x in paired_parameters]
high_ell=[x.upper() for x in paired_parameters]

def old_height_lower(ell):
    if not ell>0 or not ell<paired_parameters[-1]:return zero
    start=max(0,bisect_left(high_ell,ell.lower())-1)
    finish=min(len(paired_height),bisect_right(low_ell,ell.upper())+1)
    assert start<finish
    return min(paired_height[start:finish])

'''
    prefix=replace_line(prefix,marker,height_table+marker)
    prefix=replace_line(prefix,'cell_deficit = (own_height / value.lower()).upper()','cell_deficit = max(zero,((own_height - old_height_lower(l)) / value.lower()).upper())')
    ns={'__file__':old['__file__'],'__name__':'__main__'}
    with contextlib.redirect_stdout(io.StringIO()):exec(compile(prefix,old['__file__'],'exec'),ns)
    return ns,prefix,old
