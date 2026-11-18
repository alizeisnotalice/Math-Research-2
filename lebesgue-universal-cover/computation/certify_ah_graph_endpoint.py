"""Check the exact graph-value buffer at each endpoint before applying contact ordering."""
import inspect
from fractions import Fraction as F
from certify_ah_paired_height import geometry as base_geometry
from certify_ah_margin_envelope import replace_line


def geometry(archive,extra,*,far,negative,eps):
    original=inspect.getsource(base_geometry);marker="    ns={'__file__':old['__file__'],'__name__':'__main__'}"
    before='''    assert l>lmin+margin and l<lmax-margin,(a,l,margin)
    assert pl>lmin+margin and pl<lmax-margin,(a,pl,margin)'''
    after='''    E=max(deficit_by_bin[k] for k in used)
    for ell in (l,pl):
        assert lmin<ell<lmax, (a,ell)
        gg=graph_hull(ell)
        assert graph_value(lmin)>gg+2*E and graph_value(lmax)<gg-2*E,(a,ell,E)'''
    extra_code='    prefix=replace_line(prefix,'+repr(before)+','+repr(after)+')\n'
    original=replace_line(original,marker,extra_code+marker);scope=dict(base_geometry.__globals__);exec(compile(original,'exact_endpoint_geometry','exec'),scope);return scope['geometry'](archive,extra,far=far,negative=negative,eps=eps)
