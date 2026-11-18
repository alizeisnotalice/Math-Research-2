"""Certify a deep chord using the projected partial derivatives of both cap-point charts."""
import argparse,contextlib,inspect,io,json,sys
from fractions import Fraction as F
from pathlib import Path
from hashlib import sha256
from flint import ctx
import certify_ah_deep_single_line as single
from certify_ah_fine_pivot import make_geometry as old_geometry,fine_verify
from certify_ah_margin_envelope import replace_line
from certify_ah_adaptive_bracket import verify as original_verify


def geometry(archive,extra):
    code=inspect.getsource(old_geometry)
    marker="    ns={'__file__':str(src),'__name__':'__main__'}";code=replace_line(code,marker,'''    full=replace_line(full,'assert hmax<A("0.000012")','assert hmax<A("0.000040")')
    full=replace_line(full,'F("0.07"), F("0.43"), F("0.00065")','F("0.085"), F("0.43"), F("0.0016")')
    full=replace_line(full,'pbound = F("0.072")','pbound = F("0.089")')
'''+marker)
    nd=dict(old_geometry.__globals__);exec(compile(code,'projected_geometry','exec'),nd);return nd['make_geometry'](archive,extra)


def bracket(ns,bins):
    _,code=fine_verify(ns,[])
    code=replace_line(code,'    def projected(k,u):','''    def projected(k,u):
        if k in (2,4):
            xl,xt=jet_bounds[k]
            return min(tol[k].upper(),(A(F(2,3))*abs(dot(xl,u))+abs(dot(xt,u))).upper())''')
    marker='            v=vertices_for(a,b);';code=replace_line(code,marker,'''            ca,cb=[ns['chart'](ns['add'](ns['z'],ns['mul'](A(t),D))) for t in (a,b)]
            l,t=ca[0].union(cb[0]),ca[1].union(cb[1])
            lp=ns['coord_hull'](ns['psi'],l)
            jet_bounds={}
            for k,ell,ee in ((4,l,t),(2,lp,-t)):
                x,lc,nc=ns['chart_bounds'](ell+(-A(F(2,3))*E).union(A(F(2,3))*E),ee+(-E).union(E))
                xl,xt=ns['derivs'](x,lc,nc)
                jet_bounds[k]=(xl,xt) if k==4 else (ns['reflection'](xl),ns['reflection'](xt))
'''+marker)
    code=replace_line(code,'elif j<35: turn,sign,cells=300,1,75','elif j<34: turn,sign,cells=300,1,75\n        elif j==34: turn,sign,cells=0,1,75')
    code=replace_line(code,'(F(1,20) if j>=40 else F(1))','F(1,20)')
    stage=dict(original_verify.__globals__);exec(compile(code,'projected_brackets','exec'),stage);return stage['verify'](ns,bins),code


def green(archive,p,d):
    code=inspect.getsource(single.green);code=replace_line(code,'F("0.07")','F("0.085")');nd=dict(single.green.__globals__);exec(compile(code,'longer_green','exec'),nd);return nd['green'](archive,p,d)


def run(archive,inp,out,bins):
    extra=[tuple(map(F,r)) for r in json.loads(inp.read_bytes())]
    if bins is not None:
        ns,_=geometry(archive,extra);bracket(ns,bins);return
    code=inspect.getsource(single.run);code=replace_line(code,"0<d<=F('.0006')","F('.00083')<d<F('.002')");code=replace_line(code,"hi<F('.844093390960250720226452027837524452015')","hi<F('.844093332639707301922026011926669949130')")
    stage={**single.run.__globals__,'make_geometry':geometry,'fine_verify':bracket,'green':green,'__file__':__file__};exec(compile(code,'projected_chord','exec'),stage);stage['run'](archive,inp,out)

if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('archive',type=Path);ap.add_argument('input',type=Path);ap.add_argument('output',type=Path);ap.add_argument('--bins',nargs='+',type=int);opt=ap.parse_args();sys.path.insert(0,str(opt.archive.resolve()));run(opt.archive.resolve(),opt.input,opt.output,opt.bins)
