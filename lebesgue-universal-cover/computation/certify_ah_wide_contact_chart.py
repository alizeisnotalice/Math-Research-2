"""Extend the paired-arc window and reverify all chart inequalities."""
import inspect
from certify_ah_projected_jets import geometry as chart_geometry, green as chart_green, bracket
from certify_ah_margin_envelope import replace_line

WF='.10'
PB='.105'
EPS='.0026'

def geometry(archive, extra):
    text=inspect.getsource(chart_geometry)
    text=replace_line(text,'F("0.085"), F("0.43"), F("0.0016")',f'F("{WF}"), F("0.43"), F("{EPS}")')
    text=replace_line(text,'F("0.089")',f'F("{PB}")')
    text=replace_line(text,'A("0.000040")','A("0.000054")')
    marker="    nd=dict(old_geometry.__globals__);"
    before='''        ln=E + box(-eps,eps)
        tn=T + box(-eps,eps)
        x,lc,nc=chart_bounds(ln,tn)
        xl,xt=derivs(x,lc,nc)
        max_xl,max_xt=max(max_xl,norm(xl).upper()),max(max_xt,norm(xt).upper())
        norm_l, norm_t = ("0.24", "0.68")
        assert norm(xl) < A(norm_l) and norm(xt) < A(norm_t), (a, norm(xl), norm(xt))'''
    after='''        for lo_l,hi_l in ((-eps,0),(0,eps)):
            for lo_t,hi_t in ((-eps,0),(0,eps)):
                ln,tn=E+box(lo_l,hi_l),T+box(lo_t,hi_t)
                x,lc,nc=chart_bounds(ln,tn)
                xl,xt=derivs(x,lc,nc)
                max_xl,max_xt=max(max_xl,norm(xl).upper()),max(max_xt,norm(xt).upper())
                assert norm(xl)<A("0.24") and norm(xt)<A("0.68"), (a,norm(xl),norm(xt))'''
    text=replace_line(text,marker,
        '    code=replace_line(code,marker, "    full=replace_line(full,"+repr('+repr(before)+')+","+repr('+repr(after)+')+")\\n"+marker)\n'+marker)
    scope=dict(chart_geometry.__globals__)
    exec(compile(text,'wider_chart_geometry','exec'),scope)
    return scope['geometry'](archive,extra)


def green(archive,p,d):
    text=inspect.getsource(chart_green)
    text=replace_line(text,"'F(\"0.085\")'",f"""'F("{WF}")'""")
    scope=dict(chart_green.__globals__)
    exec(compile(text,'wider_green','exec'),scope)
    return scope['green'](archive,p,d)
