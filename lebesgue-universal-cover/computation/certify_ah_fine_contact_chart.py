"""Subdivide the contact-derivative tubes while preserving the same .24 and .68 norm bounds."""
from fractions import Fraction as F
import inspect
from certify_ah_first_hit_chain import geometry as seed
from certify_ah_two_hit_budgets import geometry as complete
from certify_ah_margin_envelope import ANCHORED_CUTS,replace_line


def fine_chart(archive,extra,*,far,negative,eps):
    old,src,_=seed(archive,[(F('.427'),F('.0048'))],far=F('.183'),negative=F('.210'),eps=F('.0051'))
    p=max(p for p,_ in extra)
    assert F('.419')<p<F('.430') and (p*5000).denominator==(far*1000).denominator==(negative*5000).denominator==1
    def change(prefix,replacement):return replace_line(src,next(row for row in src.splitlines() if row.startswith(prefix)),replacement)
    src=change('cuts=',ANCHORED_CUTS[:-1]+''.join(f', (F("{u}"), F("{v}"))' for u,v in extra)+']')
    src=change('wp, wf,',f'wp, wf, wmax, eps = F("{p}"), F("{far}"),F("0.43"), F("{eps}")')
    src=change('pbound = ',f'pbound = F("{negative}")')
    old_loop='tuple((eps*F(i,2),eps*F(i+1,2)) for i in range(-2,2))'
    assert src.count(old_loop)==2
    src=src.replace(old_loop,'tuple((eps*F(i,6),eps*F(i+1,6)) for i in range(-6,6))')
    src=replace_line(src,'for l,pl,used,a in parameter_boxes:','for (l,pl,used,a),(_,e_cell) in zip(parameter_boxes,physical_deficits,strict=True):')
    src=replace_line(src,'    E=max(deficit_by_bin[k] for k in used)','    E=e_cell')
    return old,src,None


def geometry(archive,extra,*,far,negative,eps):
    src=inspect.getsource(complete);scope={**complete.__globals__,'first_geometry':fine_chart}
    exec(compile(src,'fine_physical_root_chart','exec'),scope)
    return scope['geometry'](archive,extra,far=far,negative=negative,eps=eps)
