"""Extend only the paired graph and its slope/derivative tables; the reference pivot grid remains the same."""
import inspect
from fractions import Fraction as F
from certify_sigma_ah_full_chart import geometry as base
from certify_ah_taller_floor_chart import taller
from certify_ah_margin_envelope import replace_line


def geometry(archive,extra,*,sigma,far,negative,eps,graph_end=F('.48')):
    def widened(*args,**kw):
        ns,src,out=taller(*args,**kw)
        src=replace_line(src,'F("0.43"), F(',f'F("{graph_end}"), F(')
        src=replace_line(src,'for k in range(430):','for k in range(int(1000*wmax)):')
        src=replace_line(src,'assert norm(xl)<A("0.24") and norm(xt)<A("0.68")','assert norm(xl)<A("0.28") and norm(xt)<A("0.78")')
        assert F('.28')*F(2,3)+F('.78')<1
        return ns,src,out
    source=inspect.getsource(base);env={**base.__globals__,'taller':widened}
    exec(compile(source,'extend_sigma_partner_graph','exec'),env)
    ns,src,_=env['geometry'](archive,extra,sigma=sigma,far=far,negative=negative,eps=eps)
    ns['_contact_chain_slopes']=[(ell,(-m).lower()) for ell,m in ns['graph_slopes']]
    assert all(m>3 for ell,m in ns['_contact_chain_slopes'])
    return ns,src,None
