"""Continuous curvature dual for all threefold constant-width bodies and a uniform nonsymmetric neighborhood."""
import argparse,inspect,json,sys
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
from flint import arb,ctx
import certify_ah_tapered_release as base


def geometry(archive):
    src=inspect.getsource(base.run)
    marker='    tabs=[]'
    assert src.count(marker)==1
    src=src[:src.index(marker)]+'    return A,g,active,sprague,hexnames\n'
    scope=dict(vars(base))
    exec(compile(src,'same_tapered_host_rows','exec'),scope)
    A,g,active,sprague,hexnames=scope['run'](archive)
    L=arb.pi()/3
    def fold(normal):
        phi=arb.atan2(normal[1],normal[0])
        if phi<0:
            phi+=2*arb.pi()
        found=[(k,phi-k*L) for k in range(6) if k*L<phi<(k+1)*L]
        assert len(found)==1,(normal,phi)
        k,x=found[0]
        return (-1)**k,x
    c=next(u for name,u,h in active if name=='I-E3')
    other=next(u for name,u,h in active if name=='K-C3')
    sc,xc=fold(c)
    s5,x5=fold(other)
    assert s5==sc and (x5-xc).contains(0)
    rows=[(name,*fold(u),h) for name,u,h in active if name not in hexnames|{'I-E3','K-C3','taper_64'} and not name.startswith(sprague)]
    return A,L,sc,xc,rows


def generate(archive,out):
    import numpy as np
    from scipy.optimize import minimize
    A,L,sc,xc,rows=geometry(archive)
    L=float(L)
    v=(np.arange(2048)+.5)*L/2048
    def kernel(x):
        return np.sin(np.minimum(v,x))*np.sin(L-np.maximum(v,x))/np.sin(L)
    cap=sc*kernel(float(xc))
    moment=np.sin(v)-np.sin(L-v)
    weights={}
    for name,s,x,h in rows:
        direct=s*kernel(float(x))
        def objective(z):
            return np.mean(abs(direct-z[0]*cap-z[1]*moment))
        z=minimize(objective,[max(0.,np.mean(direct)/np.mean(cap)),0],bounds=((0,3),(-2,2)),method='Nelder-Mead',
                   options={'xatol':1e-8,'fatol':1e-10}).x
        weights[name]=[f'{z[0]:.7f}',f'{z[1]:.7f}']
    out.write_text(json.dumps(weights,indent=2)+'\n')


def prove(archive,source):
    from certify_b1_curvature_dual import positive_area
    A,L,sc,xc,rows=geometry(archive)
    weights=json.loads(source.read_bytes())
    assert set(weights)=={name for name,*_ in rows}
    sinL,cosL=L.sin_cos()
    def coefficients(x,left):
        if left:
            return (L-x).sin()/sinL,A(0)
        return -x.sin()*cosL/sinL,x.sin()
    data=[]
    for name,s,x,h in rows:
        lam,nu=map(A,weights[name])
        assert lam>=0
        ends=sorted((A(0),x,xc,L),key=lambda y:float(y.mid()))
        integral=A(0)
        for l,r in zip(ends,ends[1:]):
            assert l<r
            mid=(l+r)/2
            assert (mid<x or mid>x) and (mid<xc or mid>xc)
            ka,kb=coefficients(x,mid<x)
            ca,cb=coefficients(xc,mid<xc)
            a=s*ka-lam*sc*ca-nu*(1+cosL)
            b=s*kb-lam*sc*cb+nu*sinL
            integral+=positive_area(a,b,l,r)+positive_area(-a,-b,l,r)
        slack=(h-A('.5')-integral/2).lower()
        radius=(slack/(1+lam)).lower()
        assert slack>A('.00088') and radius>A('.00075'),(name,slack,radius)
        data.append((name,slack,radius))
    worst=min(data,key=lambda z:z[2])
    normfactor=1+2/A(3).sqrt()
    kernelmass=2/A(3).sqrt()-1
    assert normfactor*2*A('.00052')/3<A('.00075')
    assert kernelmass*2*A('.0072')/3<A('.00075')
    return dict(host='ah-tapered-release.json',all_threefold_constant_width_bodies_fit=True,
                allowed_Hausdorff_distance_to_own_threefold_average='0.00075',
                allowed_normalized_120degree_rotation_discrepancy='0.001125',
                allowed_rotation_discrepancy_about_some_center='0.00052',
                allowed_Linfty_angular_curvature_discrepancy='0.0072',signed_cap_generator=int(sc),
                active_rows=len(data),radius_limiting_row=worst[0],proven_radius_lower=str(worst[2]),
                threefold_noncap_min_slack=str(min(r[1] for r in data)),integrated_on_continuous_sector=True,
                weights_sha256=sha256(source.read_bytes()).hexdigest(),base_sha256=sha256(Path(base.__file__).read_bytes()).hexdigest(),
                verifier_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),universal_placement_proved=False)

if __name__=='__main__':
    cli=argparse.ArgumentParser()
    cli.add_argument('archive',type=Path)
    cli.add_argument('weights',type=Path)
    cli.add_argument('output',type=Path,nargs='?')
    cli.add_argument('--generate',action='store_true')
    opt=cli.parse_args()
    ar=opt.archive.resolve()
    sys.path.insert(0,str(ar))
    if opt.generate:
        generate(ar,opt.weights)
    else:
        result=prove(ar,opt.weights)
        opt.output.write_text(json.dumps(result,indent=2)+'\n')
        print(result)
