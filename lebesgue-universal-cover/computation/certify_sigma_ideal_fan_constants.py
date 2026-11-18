"""Directed full-reference proof of the sector ordering and the 1.07 support transfer at the changed Gibbs angle."""
import argparse,json,sys
from hashlib import sha256
from fractions import Fraction as F
from pathlib import Path
from certify_sigma_extended_graph import geometry


def run(base,inp,out):
    cfg=json.loads(inp.read_bytes());p,d=map(F,cfg['line']);far,negative,eps=map(F,cfg['physical_parameters']);ns,src,_=geometry(base,[(p,d)],sigma=cfg['sigma_degrees'],far=far,negative=negative,eps=eps)
    A,sub,dot,det=(ns[k] for k in ('A','sub','dot','cross'));orders=ns['orders'];low=gap=turn=side=A(1)
    for i in range(430):
        a,b=F(i,1000),F(i+1,1000);v=ns['centered_points'](a,b);edges=[sub(v[k],v[l]) for k,l in orders]
        for k,e in enumerate(edges):
            f=edges[(k+1)%10];side=min(side,(e[1] if k<5 else -e[1]).lower());turn=min(turn,det(e,f).lower());low=min(low,dot(e,f).lower())
        for k in range(5):
            for l in range(k):
                if (k,l) in orders:continue
                u=sub(v[k],v[l]);gap=min(gap,(1-dot(u,u)).lower())
    assert min(gap,turn,side)>0 and low>A('.75') and (A(2)/(1+low)).sqrt()<A('1.07');obj={'sigma_degrees':cfg['sigma_degrees'],'reference_closed_interval':['0','.43'],'parameter_cells':430,'strict_signed_vertical_lower':str(side),'successive_cross_lower':str(turn),'successive_dot_lower':str(low),'nonneighbor_squared_distance_gap_lower':str(gap),'positive_cone_sum_upper':str((A(2)/(1+low)).sqrt().upper()),'physical_geometry_sha256':sha256(src.encode()).hexdigest(),'script_sha256':sha256(Path(__file__).read_bytes()).hexdigest(),'flint':ns['flint_version'],'python':sys.version};out.write_text(json.dumps(obj,indent=2)+'\n');print(obj)


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('archive',type=Path);ap.add_argument('input',type=Path);ap.add_argument('output',type=Path);o=ap.parse_args();run(o.archive.resolve(),o.input,o.output)
