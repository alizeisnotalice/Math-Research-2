"""Certify the common containment of the two-parameter interpolation of four pinned poses."""
import argparse, json, sys
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
from flint import arb, arb_series, ctx, __version__ as flint_version

def run(comp, spec, report):
    sys.path.insert(0,str(comp));from deformed_reuleaux_support import body, environment
    from certify_tapered_forced_support import rational
    entry=json.loads(spec.read_bytes()); prior=json.loads(report.read_bytes())
    assert sha256(spec.read_bytes()).hexdigest()==prior['config_sha256']
    rows=entry['bodies']; assert len(rows)==7 and prior['family']==rows
    corners=rows[3:]
    ws=sorted(set(F(k[0]) for k in corners)); es=sorted(set(F(k[1]) for k in corners))
    assert len(ws)==len(es)==2 and {(F(z[0]),F(z[1])) for z in corners}=={(w,e) for w in ws for e in es}
    assert len({tuple(z[2:]) for z in corners})==1 and entry['translation_rule']=='pin_hexagon_faces_0_and_2'
    assert entry['family_interval']=={'first':1,'last':2}
    ns=environment(str(comp),'1.5497');ctx.prec=192
    W=ns['box'](*ws); E=ns['box'](*es); body(ns,W,E)
    def second(variable):
        w=arb_series([W,int(variable=='w')]);e=arb_series([E,int(variable=='e')])
        add,mul,reflect=ns['add'],ns['mul'],ns['reflection']
        l,t=ns['chart'](add(ns['z'],mul(w,ns['D'])))
        p=ns['psi'](l);t=t-e
        vs=[add(ns['l0'],mul(l,ns['el'])),reflect(add(ns['l0'],mul(p,ns['el']))),
            reflect(ns['X'](p,-t)),add(ns['mid'],mul(t,ns['en'])),ns['X'](l,t)]
        result=[2*(abs(x[2]).upper()**2+abs(y[2]).upper()**2).sqrt() for x,y in vs]
        assert all(z.is_finite() for z in result)
        return result
    bounds=[second(v) for v in ('w','e')]
    ds=[b-a for a,b in (ws,es)]
    # Bounds for the support acceleration after pinning are at most three
    # times the maximum acceleration of the five unrotated vertices.
    r=sum(3*max(z.upper() for z in b)*rational(d*d/8) for b,d in zip(bounds,ds))
    radius=F(14,10**11)
    assert r<rational(radius),(str(r),bounds)
    # Original six-normal container gives C* plus corners Euclidean norm <1.
    # The host certificate checks five supports <=.5 and the remaining <=.500000301
    # for the old family. Check the four new ones independently.
    assert all(arb(z) < arb('.500001') for a in prior['six_body_supports'] for z in a['six_supports'])
    # In fact this bounding hexagon lies in B(0,.58), giving perimeter <=2*pi*.58<4.
    assert arb('.500001')/(arb.pi()/6).cos() < arb('.58')
    assert 2*arb.pi()*arb('.58')<4
    upper=F(prior['unrelaxed_upper'])+4*radius+4*radius**2
    assert upper<F(prior['B'])
    return dict(rectangle=[list(map(str,v)) for v in (ws,es)],vertex_acceleration=[[str(z) for z in a] for a in bounds],
        acceleration_radius=str(r),rational_radius=str(radius),area_strict_upper=str(upper),remaining_B=str(F(prior['B'])-upper),
        corners_sha256=sha256(spec.read_bytes()).hexdigest(),receipt_sha256=sha256(report.read_bytes()).hexdigest(),
        python_flint=flint_version,precision=ctx.prec,script_sha256=sha256(Path(__file__).read_bytes()).hexdigest())
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('comp',type=Path);p.add_argument('corners',type=Path)
    p.add_argument('certificate',type=Path);p.add_argument('output',type=Path);a=p.parse_args()
    out=run(a.comp.resolve(),a.corners,a.certificate);a.output.write_text(json.dumps(out,indent=2)+'\n');print(float(F(out['area_strict_upper'])))
