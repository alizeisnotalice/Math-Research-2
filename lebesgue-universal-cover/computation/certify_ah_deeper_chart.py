"""Deep single AH chord with a reverified larger chart tube; derivative estimates depend on tube width."""
import argparse,contextlib,inspect,io,json,sys
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
from flint import ctx
import certify_ah_deep_single_line as prior
from certify_ah_fine_pivot import fine_verify, make_geometry as initial_geometry
from certify_ah_margin_envelope import replace_line


def geometry(archive,extra):
    code=inspect.getsource(initial_geometry)
    marker="    ns={'__file__':str(src),'__name__':'__main__'}"
    # All physical caps, inverse maps, derivative bounds, trigger and graph-end conditions are recomputed.
    code=replace_line(code,marker,'''    full=replace_line(full,'assert hmax<A("0.000012")','assert hmax<A("0.000030")')
    full=replace_line(full,'F("0.43"), F("0.00065")','F("0.43"), F("0.0011")')
'''+marker)
    ns=dict(initial_geometry.__globals__);exec(compile(code,'larger_chart','exec'),ns)
    return ns['make_geometry'](archive,extra)


def run(archive,inp,out,bins=None):
    extra=[tuple(map(F,r)) for r in json.loads(inp.read_bytes())];(p,d),=extra
    if bins is not None:
        chart,full=geometry(archive,extra);fine_verify(chart,bins);return
    code=inspect.getsource(prior.run);code=replace_line(code,"0<d<=F('.0006')","F('.0006')<d<F('.002')")
    code=replace_line(code,"hi<F('.844093390960250720226452027837524452015')","hi<F('.844093374562671835178166470806391361029')")
    ns={**prior.run.__globals__,'make_geometry':geometry,'__file__':__file__};exec(compile(code,'deep_chart_verify','exec'),ns);ns['run'](archive,inp,out)

if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('archive',type=Path);ap.add_argument('input',type=Path);ap.add_argument('output',type=Path);ap.add_argument('--bins',nargs='+',type=int);opt=ap.parse_args();sys.path.insert(0,str(opt.archive.resolve()));run(opt.archive.resolve(),opt.input,opt.output,opt.bins)
