"""A physical AH chart valid at recomputed Gibbs angles, with tighter derivative subdivision for deeper cuts."""
from fractions import Fraction as F
import inspect
from certify_ah_taller_floor_chart import taller
from certify_ah_two_hit_budgets import geometry as assemble
from certify_ah_margin_envelope import replace_line


def geometry(archive,extra,*,sigma,far,negative,eps):
    def initial(archive,extra,**kw):
        ns,src,_=taller(archive,extra,**kw);src=replace_line(src,'sigma_degrees = "1.5497"',f'sigma_degrees = "{sigma}"')
        src=replace_line(src,'assert hmax<A("0.00024")','assert hmax<A("0.00035")')
        old='tuple((eps*F(i,6),eps*F(i+1,6)) for i in range(-6,6))'
        assert src.count(old)==2;src=src.replace(old,'tuple((eps*F(i,10),eps*F(i+1,10)) for i in range(-10,10))')
        return ns,src,None
    env={**assemble.__globals__,'first_geometry':initial}
    exec(compile(inspect.getsource(assemble),'sigma_large_cap_specific_geometry','exec'),env)
    return env['geometry'](archive,extra,far=far,negative=negative,eps=eps)
