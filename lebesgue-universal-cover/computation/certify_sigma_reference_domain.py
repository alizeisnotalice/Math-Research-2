"""Certify that extending the partner derivative graph leaves every actual first hit inside the reference grid."""
import argparse,json,sys
from fractions import Fraction as F
from pathlib import Path
from hashlib import sha256
from certify_sigma_extended_graph import geometry


def run(archive,inp,out):
    cfg=json.loads(inp.read_bytes());p,d=map(F,cfg['line']);far,negative,eps=map(F,cfg['physical_parameters']);ns,src,_=geometry(archive,[(p,d)],sigma=cfg['sigma_degrees'],far=far,negative=negative,eps=eps)
    A=ns['A'];r0=ns['chart'](ns['add'](ns['z'],ns['mul'](A(ns['knots'][-1]),ns['D'])))[0];upper=ns['coord_hull'](ns['psi'],r0)
    gap=min(min((point-r0).lower(),(upper-point).lower()) for ell,_ in ns['physical_deficits'] for point in (ell,ns['coord_hull'](ns['psi'],ell)))
    assert gap>0
    out.write_text(json.dumps({'line':cfg['line'],'sigma':cfg['sigma_degrees'],'reference_grid_closes_at':str(ns['knots'][-1]),'lower_chart_coordinate':str(r0),'upper_chart_coordinate':str(upper),'min_signed_gap_to_either_end':str(gap),'physical_rows':len(ns['physical_deficits']),'geometry_sha256':sha256(src.encode()).hexdigest(),'parameters_sha256':sha256(inp.read_bytes()).hexdigest(),'flint':ns['flint_version'],'python':sys.version},indent=2)+'\n');print(gap)


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('archive',type=Path);ap.add_argument('input',type=Path);ap.add_argument('output',type=Path);o=ap.parse_args();run(o.archive.resolve(),o.input,o.output)
