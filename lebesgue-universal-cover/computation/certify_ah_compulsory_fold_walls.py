"""Force the paired AH boundary down to the ideal compatibility fold and evaluate both sharp endpoint walls."""
import argparse,json,sys
from fractions import Fraction as F
from pathlib import Path
from hashlib import sha256
from certify_ah_two_hit_budgets import geometry


def run(archive,fold,obstruction,out):
    a=json.loads(fold.read_bytes());b=json.loads(obstruction.read_bytes());assert b['fold_certificate_sha256']==sha256(fold.read_bytes()).hexdigest()
    ns,_,_=geometry(archive,[(F('.427'),F('.0048'))],far=F('.183'),negative=F('.210'),eps=F('.0051'))
    A=ns['A'];dot,sub,mul,add=map(ns.get,('dot','sub','mul','add'));l,r=map(F,a['unique_fold_exact_enclosure']);w=ns['box'](l,r)
    d=ns['direct'];n=ns['normal'];z=ns['z'];L=ns['length'];vertices=[(A(0),1/ns['r3']),(A('.5'),A('.5')/ns['r3']),(A('.5'),-A('.5')/ns['r3']),(A(0),-1/ns['r3']),(-A('.5'),-A('.5')/ns['r3']),(-A('.5'),A('.5')/ns['r3'])]
    violation=[(-dot(n,sub(x,z))).lower() for x in vertices];best=max(violation);k=violation.index(best);assert best>0
    point=ns['reflection'](ns['centered_points'](l,r)[2]);lam=dot(d,sub(point,z))/L;eta=dot(n,sub(point,z))/L
    # The original paired graph has increasing projected coordinate in this extension.
    low=A(1)
    for i in range(50):
        x=F('.428')+F(i,25000);y=x+F(1,25000);v=ns['reflection'](ns['point_jets'](ns['box'](x,y))[1][2])
        low=min(low,(dot(d,v)/L).lower())
    assert low>0 and lam>0 and eta>0
    tail=json.loads((archive/'ah-compulsory-tail.json').read_bytes())
    obj={'w_fold':str(w),'paired_projection_lambda':str(lam),'paired_height_eta':str(eta),'paired_projection_derivative_lower':str(low),'hex_vertex_for_straight_escape':k,'hex_vertex_coordinates':[str(x) for x in vertices[k]],'vertex_excess_over_straight_support':str(best),'previous_tail_sha256':sha256((archive/'ah-compulsory-tail.json').read_bytes()).hexdigest(),'previous_tail_receipt_keys':sorted(tail.keys()),'fold_certificate_sha256':sha256(fold.read_bytes()).hexdigest(),'obstruction_certificate_sha256':sha256(obstruction.read_bytes()).hexdigest(),'source_sha256':sha256(Path(__file__).read_bytes()).hexdigest(),'python':sys.version,'flint':ns['flint_version'],'bits':128}
    out.write_text(json.dumps(obj,indent=2)+'\n');print(str(lam),str(eta),k,str(best))


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('archive',type=Path);ap.add_argument('fold',type=Path);ap.add_argument('obstruction',type=Path);ap.add_argument('output',type=Path);opt=ap.parse_args();run(opt.archive.resolve(),opt.fold,opt.obstruction,opt.output)
