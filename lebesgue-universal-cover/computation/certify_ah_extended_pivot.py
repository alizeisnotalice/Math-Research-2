"""Verify a longer AH pivot by adaptive angular bracketing and true-boundary area."""
import argparse
import contextlib
from fractions import Fraction as F
from hashlib import sha256
import io
import json
from pathlib import Path
import sys

from flint import arb, ctx
from certify_ah_adaptive_bracket import verify
from certify_ah_margin_envelope import ANCHORED_CUTS, coefficients, integral, replace_line, upper_envelope
from certify_ah_three_line_cut import dyadic, fit_archived_witnesses, format_area, specialized_source


def make_geometry(archive, extra):
    path=archive/'certify_ah_two_line_cap.py'
    full=specialized_source(path,F('.404'),F('.00046'))
    cuts=ANCHORED_CUTS[:-1]+''.join(f', (F("{p}"), F("{d}"))' for p,d in extra)+']'
    full=replace_line(full,ANCHORED_CUTS,cuts)
    top=max(p for p,_ in extra); assert F('.410')<top<F('.43') and top*1000==int(top*1000)
    full=replace_line(full,'wp, wf, wmax, eps = F("0.410")',f'wp, wf, wmax, eps = F("{top}")')
    ns={'__file__':str(path),'__name__':'__main__'}
    with contextlib.redirect_stdout(io.StringIO()):
        exec(compile(full.split('angular_cells=32')[0],str(path),'exec'),ns)
    return ns,full


def area(ar):
    import audit_optimized_eh_chord as g
    from audit_eh_continuous_placement import polygon
    from certify_eh_polyline_obstruction import clip_geometry
    z,D,n,L,A=(ar[k] for k in ('z','D','normal','length','A'))
    lines=coefficients(ar); pmax=max(p for p,_ in ar['cuts']); end=A(pmax+F('.001'))
    assert ar['wp']==pmax and end<1
    zero=(arb(0),arb(0),-1)
    for b,a,_ in lines[3:]:
        assert lines[0][0]>b>0 and lines[0][1]>a
    before=upper_envelope(lines[:3]+[zero],end); after=upper_envelope(lines+[zero],end)
    assert after[0][2]==before[0][2]==0 and after[-1][2]==before[-1][2]==-1
    lo,hi=after[1][3],after[-1][3]
    knots=[r[3] for r in after[1:]]
    for *_,w in before[1:]:
        assert w<lo or w>hi or lo<w<hi
        if lo<w<hi: knots.append(w)
    g.SIGMA=A('1.5497')*g.PI/180; inner,_=polygon('0.0005522',32)
    dirs,offsets,*_=clip_geometry(['0.000555','0.0005522','0.0075'])
    inner.extend((f'R126_EH_{i}',u,h) for i,(u,h) in enumerate(zip(dirs,offsets)))
    # The straight Ta-Z side is exactly y >= 0. Both envelopes include y=0 as a supporting line.
    assert sum(name=='Ta-Z' for name,_,_ in inner)==1
    points=[]
    for w in knots:
        assert 0<w<1
        for subgroup in (lines[:3],lines):
            ys=[a-b*w for b,a,_ in subgroup+[zero]]
            y=max(t.lower() for t in ys).union(max(t.upper() for t in ys))
            points.append(g.plus(g.plus(z,g.scale(w,D)),g.scale(L*y,n)))
    clear=min((c-g.dot(u,x)).lower() for name,u,c in inner if name!='Ta-Z' for x in points)
    assert clear>0
    gain=L**2*(integral(after,end)-integral(before,end)); assert gain>0
    return {'old_active':[r[2] for r in before],'new_active':[r[2] for r in after],'new_breakpoints':[str(r[3]) for r in after[1:]],'other_curve_clearance_lower':str(clear),'area_gain_exact':[str(dyadic(q)) for q in (gain.lower(),gain.upper())]},gain,inner


def run(path,inp,out,bins):
    extra=[tuple(map(F,r)) for r in json.loads(inp.read_bytes())]
    ar,full=make_geometry(path,extra); assert all(0<p<F('.43') and 0<d<F('.00060') for p,d in extra)
    checked=verify(ar,bins)
    if bins is not None: return
    ctx.prec=160; shape,gain,inner=area(ar)
    prior_path=path/'ah-three-line-cut-certificate.json'; prior=json.loads(prior_path.read_bytes()); a,b=map(F,prior['area_exact_interval']); lo,hi=a-dyadic(gain.upper()),b-dyadic(gain.lower())
    assert lo<hi<a and hi<F('.844093393709073433643894640415099725075')
    out.write_text(json.dumps({'lines':json.loads(inp.read_bytes()),'derived_source_sha256':sha256(full.encode()).hexdigest(),'R167_area_sha256':sha256(prior_path.read_bytes()).hexdigest(),'input_sha256':sha256(inp.read_bytes()).hexdigest(),'script_sha256':sha256(Path(__file__).read_bytes()).hexdigest(),'verification_module_sha256':sha256(Path(__file__).with_name('certify_ah_adaptive_bracket.py').read_bytes()).hexdigest(),'python':sys.version,'flint':ar['flint_version'],'geometry_and_placement_bits':128,'area_bits':160,'bins':checked,'shape':shape,'area_exact_interval':[str(lo),str(hi)],'area_decimal_interval':format_area(lo,hi),'R150_R156_witnesses_at_identity':fit_archived_witnesses(path,ar,inner)},indent=2)+'\n')

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__); parser.add_argument('archive',type=Path); parser.add_argument('lines',type=Path); parser.add_argument('output',type=Path); parser.add_argument('--bins',nargs='+',type=int)
    opt=parser.parse_args(); sys.path.insert(0,str(opt.archive.resolve())); run(opt.archive.resolve(),opt.lines,opt.output,opt.bins)
