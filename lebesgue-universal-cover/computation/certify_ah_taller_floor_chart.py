"""Check a taller physical cap directly using the true lower paired arc and existing chart gates."""
import inspect
from certify_ah_fine_floor_chart import floor_chart
from certify_ah_two_hit_budgets import geometry as assemble
from certify_ah_margin_envelope import replace_line


def taller(archive,extra,*,far,negative,eps):
    ns,src,_=floor_chart(archive,extra,far=far,negative=negative,eps=eps)
    src=replace_line(src,'assert hmax<A("0.00016")','assert hmax<A("0.00024")')
    return ns,src,None


def geometry(archive,extra,*,far,negative,eps):
    env={**assemble.__globals__,'first_geometry':taller}
    exec(compile(inspect.getsource(assemble),'taller_physical_domain','exec'),env)
    return env['geometry'](archive,extra,far=far,negative=negative,eps=eps)
