"""Certify a genuine constant-width obstruction to every AH pivot beyond the ideal compatibility fold."""
import argparse,sys,json
from fractions import Fraction as F
from flint import arb
from pathlib import Path
from hashlib import sha256
from certify_ah_two_hit_budgets import geometry


def run(archive,fold,out):
    bound=json.loads(fold.read_bytes());base=archive/'certify_ah_ideal_root_fold.py'
    assert bound['source_sha256']==sha256(base.read_bytes()).hexdigest()
    wl,wh=map(F,bound['unique_fold_exact_enclosure'])
    ns,_,_=geometry(archive,[(F('.427'),F('.0048'))],far=F('.183'),negative=F('.210'),eps=F('.0051'))
    A=ns['A'];R=1/ns['r3'];V=(A(0),R);dot,sub,mul,perp,norm=map(ns.get,('dot','sub','mul','perp','norm'))
    rot,trig,sup,reflect=map(ns.get,('rotate','trig','supports','reflection'));hp=list(map(ns['av'],ns['d']['n']));mp=list(map(ns['av'],ns['d']['m']));half=A('.5')
    other=[(k,s) for k in range(0,360,60) for s in (1,-1) if (k,s) not in [(0,1),(240,-1)]]
    def data(a,b):
        v=ns['centered_points'](a,b);e=[sub(v[i],v[j]) for i,j in ns['orders']]
        q=reflect(v[2]);der=reflect(ns['point_jets'](ns['box'](a,b))[1][2]);qn=mul(-1/norm(der),perp(der))
        assert der[0]<0 and der[1]<0
        projected=dot(qn,sub(V,q)).lower()
        misses={}
        for k,s in other:
            violations=[(sup(v,e,trig(A(-k)),(mp[j][0],s*mp[j][1]))-half).lower() for j in (3,5)]
            best=max(violations);assert best>0,(a,k,s,violations);misses[f'{k},{s}']=(best,3 if violations[0]>violations[1] else 5)
        assert projected>0
        return v,projected,misses
    minproj=mine=A(1);allrows=[]
    for i in range(50):
        a=F('.428')+F(i,25000);b=a+F(1,25000);v,proj,misses=data(a,b);minproj=min(minproj,proj);mine=min(mine,*(r[0] for r in misses.values()))
        allrows.append({'w':[str(a),str(b)],'chosen_cap':{k:r[1] for k,r in misses.items()},'worst_positive_cap_violation':str(min(r[0] for r in misses.values())),'reflected_triangle_minus_true_tangent':str(proj)})
    # A concrete nonuniversal host.  The same proof lets the perturbation shrink with p-w.
    w,p,d,t=map(F,('.429','.430','.006','.00000001'));assert wh<w<p
    v,proj,misses=data(w,w)
    direction=ns['add'](mul(A(d).cos(),ns['direct']),mul(A(d).sin(),ns['normal']));cut=mul(-1,perp(direction));pivot=ns['add'](ns['z'],mul(A(p),ns['D']))
    identity=(dot(cut,sub(v[4],pivot))-2*R*A(t)).lower()
    reflected=(A(t)*proj).lower()
    others=(min(r[0] for r in misses.values())-2*R*A(t)).lower()
    coeff=[sum((-1)**k*dot(v[i],(perp if tangent else lambda x:x)(hp[k])) for k,i in enumerate((2,3,0))) for tangent in (False,True)]
    peak=(sum(x*x for x in coeff).sqrt()-half).upper()
    away=min(-peak,arb(bound['strict_negativity_margin_in_cells']).lower())-6*R*A(t)
    sw=sum((-1)**k*dot(v[i],perp(hp[k])) for k,i in enumerate((1,3,0)))
    end_slope=-((1-A(t))*sw+3*R*A(t)+A('1.5')*(A('.25')*ns['angle_unit']).tan()).upper()
    assert min(identity,reflected,others,away,end_slope)>0
    doc={'fold_certificate_sha256':sha256(fold.read_bytes()).hexdigest(),'uniform_parameter_window':['.428','.430'],'min_other_dihedral_violation_of_actual_m3_or_m5':str(mine),'min_reflected_triangle_vertex_above_true_paired_tangent':str(minproj),'uniform_cells':allrows,'explicit_witness':{'reference_w':str(w),'pivot_w':str(p),'AH_angle_radians':str(d),'Reuleaux_triangle_Minkowski_weight':str(t),'identity_cut_violation_lower':str(identity),'QA_paired_tangent_violation_lower':str(reflected),'other_ten_hexagonal_cap_violations_lower':str(others),'interior_root_exclusion_margin_away_from_0_and_60':str(away),'endpoint_sign_slope_margin':str(end_slope)},'python':sys.version,'flint':ns['flint_version'],'bits':128,'source_sha256':sha256(Path(__file__).read_bytes()).hexdigest()}
    out.write_text(json.dumps(doc,indent=2)+'\n');print(str(minproj),str(mine),doc['explicit_witness'])


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('archive',type=Path);ap.add_argument('fold',type=Path);ap.add_argument('output',type=Path);opt=ap.parse_args();run(opt.archive.resolve(),opt.fold,opt.output)
