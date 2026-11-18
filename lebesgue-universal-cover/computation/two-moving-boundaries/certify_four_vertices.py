"""Check the inherited forced arcs over every angle slab and the quantitative four-vertex consequence."""
from argparse import ArgumentParser
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import json,sys
from flint import arb,ctx
HERE=Path(__file__).resolve().parent

def run(comp,previous):
    sys.path.insert(0,str(comp))
    from certify_sigma_mandatory_arcs import anchors,certify_cell,g,degree
    from certify_three_sigma_area_budget import pack,exact
    ctx.prec=192
    prior=json.loads(previous.read_bytes())
    assert prior['complete'] and prior['interval']==['3/2','77/50'] and prior['panels']==200
    assert prior['closed_root_boxes']==730877 and prior['closed_TI_boxes']==16774
    rows=[]
    for i,r in enumerate(prior['panels_data']):
        a,z=map(F,r['degrees'])
        assert a==F('1.50')+F(i,5000) and z==a+F(1,5000)
        record={}
        for label in ('GK','QY'):
            n,stats=certify_cell(label,a,z)
            record[label]={'angle_boxes':n,'positive_lower':pack(min(v for k,v in stats.items() if 'negative' not in k).lower())}
            assert min(v for k,v in stats.items() if 'negative' not in k)>0
        # The C anchor of QY is the exact same intersection n4 and m3 as E3.
        # Source definitions give equality, and both evaluations must overlap.
        pts=anchors('QY')[0];n,m,*_=g.setup();e3=g.line_pair(n[4],m[3])
        assert all(a.overlaps(b) for a,b in zip(e3,pts['C'],strict=True))
        def proj(point,vec):
            return g.minus(point,g.scale(g.dot(point,vec)-g.HALF,vec))
        stored={s:g.line_pair(n[i],n[j]) for s,i,j in [('B1',0,1),('F1',2,3),('D1',4,5)]}
        stored['C3']=g.line_pair(n[0],m[5]);stored['E3']=e3
        f3=g.line_pair(n[3],m[2]);gg=g.circle_line(f3,n[5],1)
        stored.update(K=proj(f3,m[5]),P=proj(gg,n[2]),Q=proj(e3,n[1]),G=gg,E=anchors('TI')[0]['E'])
        triples=[('C3','B1','Q'),('P','F1','E'),('E3','D1','G'),('K','C3','B1')]
        coeffs={b:g.norm(g.minus(stored[a],stored[c]))/2 for a,b,c in triples}
        cross_min=min(g.cross(stored[a],stored[b]).lower() for group in triples for a,b in zip(group,group[1:]))
        assert cross_min>exact('.12') and all(g.norm(v)<1 for v in stored.values())
        record['coeff_upper']={k:pack(v.upper()) for k,v in coeffs.items()}
        record['radial_cross_lower']=pack(cross_min)
        rows.append(record)
    gap=F(prior['minimum_minus_B']['exact'][0])
    delta=F('0.00000009387515')
    loss=2*arb(3).sqrt()*exact(delta)+arb.pi()*exact(delta)**2
    assert loss.upper()<exact(gap).lower()
    name4=('B1','F1','D1','C3')
    coeff={k:max(F(r['coeff_upper'][k]['exact'][1]) for r in rows) for k in name4}
    d2=F('.0000002306')
    loss4=sum(coeff.values())*d2+d2*d2/2
    assert loss4<gap
    files=('certify_sigma_mandatory_arcs.py','certify_sigma_paired_rigidity.py','certify_gibbs_area_refined.py')
    return dict(complete=True,four_vertex_displacement=dict(delta=str(d2),linear_coefficients_upper={k:pack((v,v),22) for k,v in coeff.items()},
                area_loss_upper=str(loss4),positive_rational_remainder=str(gap-loss4)),constrained_vertices={'B1':'n0 & n1','F1':'n2 & n3','D1':'n4 & n5','C3':'n0 & m5'},
                mandatory_anchor='E3=n4 & m3 is QY:C',interval=prior['interval'],gap_lower=pack((gap,gap),36),
                delta=str(delta),Steiner_area_maximum=pack(loss,38),residual_lower=pack((exact(gap)-loss).lower(),38),
                angular_coverage={k:sum(r[k]['angle_boxes'] for r in rows) for k in ('GK','QY')},
                slabs=rows,continuous_certificate_sha256=sha256(previous.read_bytes()).hexdigest(),
                sources={n:sha256((comp/n).read_bytes()).hexdigest() for n in files},
                verifier_sha256=sha256(Path(__file__).read_bytes()).hexdigest())

if __name__=='__main__':
    p=ArgumentParser()
    for arg in ('comp','prior','out'):p.add_argument(arg,type=Path)
    a=p.parse_args();r=run(a.comp.resolve(),a.prior.resolve());a.out.write_text(json.dumps(r,indent=1)+'\n')
    print(r['angular_coverage'],r['delta'],r['residual_lower']['decimal'],flush=True)
