"""Refine the reference parameter bins near the small-angle fold and certify a longer pivot."""
from fractions import Fraction as F
import argparse, contextlib, io, json, sys
from hashlib import sha256
from pathlib import Path
from flint import ctx
from certify_ah_adaptive_bracket import verify as original_verify
from certify_ah_extended_pivot import area
from certify_ah_margin_envelope import ANCHORED_CUTS, replace_line
from certify_ah_three_line_cut import dyadic,fit_archived_witnesses,format_area,specialized_source


def make_geometry(path,extra):
    src=path/'certify_ah_two_line_cap.py'; full=specialized_source(src,F('.404'),F('.00046'))
    cuts=ANCHORED_CUTS[:-1]+''.join(f', (F("{p}"), F("{d}"))' for p,d in extra)+']'
    full=replace_line(full,ANCHORED_CUTS,cuts)
    top=max(p for p,_ in extra); assert F('.419')<top<F('.430') and top*5000==int(top*5000)
    full=replace_line(full,'wp, wf, wmax, eps = F("0.410")',f'wp, wf, wmax, eps = F("{top}")')
    full=replace_line(full,'ref_l=[chart(add(z,mul(A(F(i,100)),D)))[0] for i in range(int(100*wmax)+1)]',
                      'knots=[F(i,100) for i in range(41)]+[F(i,1000) for i in range(401,431)]\nref_l=[chart(add(z,mul(A(t),D)))[0] for t in knots]')
    full=replace_line(full,'for j in range(int(1000*(wp+pbound))):\n    a,b=-pbound+F(j,1000),-pbound+F(j+1,1000)',
                      'for j in range(int(5000*(wp+pbound))):\n    a,b=-pbound+F(j,5000),-pbound+F(j+1,5000)')
    ns={'__file__':str(src),'__name__':'__main__'}
    with contextlib.redirect_stdout(io.StringIO()):exec(compile(full.split('angular_cells=32')[0],str(src),'exec'),ns)
    return ns,full


def fine_verify(ar,bins):
    import inspect
    code=inspect.getsource(original_verify)
    code=replace_line(code,'seeds=high if turn else low','seeds=high if turn else low[:40]+[low[q//10]+(low[q//10+1]-low[q//10])*F(q%10,10) for q in range(400,420)]+[5231,5280,5340,5400,5450,5530,5640,5740,5900,6100,6200]\n        if j>=40: cells=20')
    code=replace_line(code,'a,b=(F(j,100)+F(i,100*cells), F(j,100)+F(i+1,100*cells))','x,y=ns["knots"][j:j+2]\n            a,b=x+(y-x)*F(i,cells),x+(y-x)*F(i+1,cells)')
    code=replace_line(code,'if widths else F(1)','if widths else (F(1,20) if j>=40 else F(1))')
    code=replace_line(code,'    def error(edges,u):', '    def projected(k,u):\n        tangent = {0:ns["el"],1:ns["reflection"](ns["el"]),3:ns["en"]}.get(k)\n        return tol[k] if tangent is None else min(tol[k].upper(),(tol[k]*abs(dot(tangent,u))).upper())\n    def error(edges,u):')
    code=replace_line(code,'(left*tol[c1]+right*tol[c2]).upper()', '(max(zero,left.upper())*projected(c1,edges[j])+max(zero,right.upper())*projected(c2,edges[nxt])).upper()')
    code=replace_line(code,'intervals.append((-upper,tol[k].upper()) if arc else (-tol[k].upper(),upper))','lower=projected(k,u).upper()\n            intervals.append((-upper,lower) if arc else (-lower,upper))')
    ns=dict(original_verify.__globals__); exec(compile(code,'fine_root_verify','exec'),ns)
    return ns['verify'](ar,bins),code


def run(path,inp,out,bins):
    extra=[tuple(map(F,l)) for l in json.loads(inp.read_bytes())]; assert all(0<p<F('.43') and 0<d<F('.00060') for p,d in extra)
    ar,full=make_geometry(path,extra); checked,source=fine_verify(ar,bins)
    if bins is not None:return
    ctx.prec=160; shape,gain,inner=area(ar)
    prior=path/'ah-three-line-cut-certificate.json'; a,b=map(F,json.loads(prior.read_bytes())['area_exact_interval']); lo,hi=a-dyadic(gain.upper()),b-dyadic(gain.lower()); assert lo<hi<F('.844093392265747203661418896791805278822')
    obj={'lines':json.loads(inp.read_bytes()),'geometry_sha256':sha256(full.encode()).hexdigest(),'fine_verifier_sha256':sha256(source.encode()).hexdigest(),'R167_area_sha256':sha256(prior.read_bytes()).hexdigest(),'input_sha256':sha256(inp.read_bytes()).hexdigest(),'script_sha256':sha256(Path(__file__).read_bytes()).hexdigest(),'python':sys.version,'flint':ar['flint_version'],'geometry_and_placement_bits':128,'area_bits':160,'bins':checked,'shape':shape,'area_exact_interval':[str(lo),str(hi)],'area_decimal_interval':format_area(lo,hi),'R150_R156_witnesses_at_identity':fit_archived_witnesses(path,ar,inner)}
    out.write_text(json.dumps(obj,indent=2)+'\n')

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__); parser.add_argument('archive',type=Path); parser.add_argument('lines',type=Path); parser.add_argument('output',type=Path); parser.add_argument('--bins',nargs='+',type=int)
    opt=parser.parse_args(); sys.path.insert(0,str(opt.archive.resolve())); run(opt.archive.resolve(),opt.lines,opt.output,opt.bins)
