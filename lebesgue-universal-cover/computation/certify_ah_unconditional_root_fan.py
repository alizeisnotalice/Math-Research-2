"""Use the free multiple of the compatibility identity on an unconditioned root bracket."""
import inspect
from certify_ah_affine_root_class import validate as base
from certify_ah_same_fan_cancellation import support_bound
from certify_ah_pinned_normals import pinned_equation,pinned_slope
from certify_ah_margin_envelope import replace_line


def validate(ns,j,left,right,turn,sign,orientation,angular=12,depth=4):
    src=inspect.getsource(base)
    src=replace_line(src,'    assert orientation==1 and condition_sign in (-1,1)',
                     '    assert orientation in (-1,1)')
    src=replace_line(src,"    anchor_terms=[(condition_sign*c,u) for c,u in root_terms(A(anchor))];host=[]",
                     '    host=[]')
    start=src.index('        f_anchor=condition_sign*')
    stop=src.index('        def visit(',start)
    src=src[:start]+'''        slope=pinned_slope(ns,v,E)
        def lower(model,minus_terms,target=zero):
            return (model-E*bound(minus_terms)).lower()
        for phi,s in ((left,orientation),(right,-orientation)):
            terms=root_terms(A(phi))
            slack=max(lower(s*ns['F_at'](v,e,phi),[(-s*c,u) for c,u in terms]),
                      (s*pinned_equation(ns,slope,A(phi))).lower())
            if slack<rm:rm,rid=slack,[i,str(phi)]
'''+src[stop:]
    src=replace_line(src,'            if gap>0:', '''            pin=pinned_equation(ns,slope,interval)
            gap=max(gap,pin.lower(),(-pin).lower())
            if gap>0:''')
    env={**base.__globals__,'support_bound':support_bound,
         'pinned_equation':pinned_equation,'pinned_slope':pinned_slope}
    exec(compile(src,'ordinary_free_root_fan','exec'),env)
    return env['validate'](ns,j,left,right,turn,sign,orientation,angular=angular,depth=depth)
