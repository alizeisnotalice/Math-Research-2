"""Classify both roots at the moving low-angle fold on a closed w neighborhood."""
from fractions import Fraction as F
from hashlib import sha256
import inspect
import json
from pathlib import Path
import sys
from flint import arb_series

import certify_two_sigma_ideal_placement as b

g=b.g
LOW,HIGH=F('.42948'),F('.42956')
LEFT,RIGHT=F('53.2'),F('54.6')
FOLD_LOW,FOLD_HIGH=F('.4295222514414'),F('.4295222514415')
CYCLE=(0,2,7,4,6)
ACTIVE=((4,True),(2,False),(6,True))


def verify():
    sigma=F('1.54')
    D=b.xy(sigma)[5]
    def coeff(w):
        p=b.anchors(D,w)
        dirs=[b.directions[i] for i in (0,2,4)]
        return [sum(b.dot(fn(d),p[i]) for d,(i,_) in zip(dirs,ACTIVE,strict=True))
                for fn in (lambda d:d,b.rot)]
    def discr(w):
        A,B=coeff(w)
        return A*A+B*B-g.HALF*g.HALF
    interval=b.ea(LOW,HIGH)
    low,high=(discr(g.exact(q)) for q in (FOLD_LOW,FOLD_HIGH))
    assert low>0 and high<0
    factors=[]
    cone=[]
    for j in range(32):
        wa,wb=[LOW+(HIGH-LOW)*F(i,32) for i in (j,j+1)]
        w=b.ea(wa,wb)
        p=b.anchors(D,w)
        midangle=b.degree(LEFT).union(b.degree(RIGHT))
        for k,(i,arc) in enumerate(ACTIVE):
            j0=CYCLE.index(i)
            v,v_next=b.sub(p[i],p[CYCLE[j0-1]]),b.sub(p[i],p[CYCLE[(j0+1)%len(CYCLE)]])
            if b.cross(v,v_next)<0:
                v,v_next=v_next,v
            assert b.cross(v,v_next)>0
            d=g.unit(midangle+k*2*g.PI/3)
            if arc:
                d=b.neg(d)
            cone += [b.cross(v,d).lower(),b.cross(d,v_next).lower()]
        a1,b1=coeff(w)
        phi=g.arb.atan2(-b1,-a1)
        assert b.degree(LEFT)<phi<b.degree(RIGHT)
        factors.append(-discr(arb_series([w,1],2))[1])
    assert min(factors)>0 and min(cone)>0
    # Reuse the full root interval guard without imposing a sign inside the
    # chamber. This also checks slope in both end neighborhoods and old caps.
    src=inspect.getsource(b.rig)
    target='scan(endcap,F(60)-endcap)'
    assert src.count(target)==1
    code=src.replace(target,'scan(endcap,F("53.2"));scan(F("54.6"),F(60)-endcap)')
    namespace=b.rig.__globals__.copy()
    exec(compile(code,'rig_outside_low_fold','exec'),namespace)
    outer=namespace['rig']
    num=0; outer_gaps=[]
    for j in range(8):
        a,c=[LOW+(HIGH-LOW)*F(k,8) for k in (j,j+1)]
        nn,bounds=outer(sigma,sigma,a,c,endcap=F('.1'),min_cell=F(1,1000000))
        num+=nn
        outer_gaps.append(bounds[-1])
    return {'complete':True,'sigma':'1.54','w_interval':list(map(str,(LOW,HIGH))),
            'chamber_degrees':list(map(str,(LEFT,RIGHT))),'discriminant_flip_interval':[str(FOLD_LOW),str(FOLD_HIGH)],
            'discriminant_at_two_bracket_ends':[str(low),str(high)],
            'affine_harmonic_formula':'F(phi)=a cos phi+b sin phi+1/2',
            'active_vertices_inverse_tip':ACTIVE,
            'minus_discriminant_w_derivative_lower':str(min(factors)),
            'angular_cone_margin_lower':str(min(cone)),
            'outside_chamber_positive_minimum':str(min(outer_gaps)),
            'outside_chamber_cells':num,
            'source_sha256':sha256(Path(__file__).read_bytes()).hexdigest()}

if __name__=='__main__':
    Path(sys.argv[1]).write_text(json.dumps(verify(),indent=2)+'\n')
