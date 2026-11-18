"""Directed interval composite-midpoint certificate for a real-space box."""
import argparse,json
from pathlib import Path
from certify import Box, D, LO, HI, ZERO, ONE, TWO, E, gaussian, gaussian_box, weight


JUMP1, JUMP2 = D('0.5'), D(1)
SHIFT1, SHIFT2 = D('0.375'), D('0.875')


def unsign(box):
    return max(box.a.copy_abs(),box.b.copy_abs())


def squared(box):
    a,b=box.a,box.b
    low=D(0) if a<=0<=b else min(LO.multiply(a,a), LO.multiply(b,b))
    return Box(low,max(HI.multiply(a,a),HI.multiply(b,b)))


def weight_cell(wlo,whi,e_lo,e_hi,crit,peak):
    if e_hi.b<=crit.a or e_lo.a>=crit.b:
        return Box(min(wlo.a,whi.a),max(wlo.b,whi.b))
    return Box(min(wlo.a,whi.a),peak.b)


def gaussian_term(lo,hi,mid,xlo,xhi,w,w_mid,first_log_weight,second_log_weight):
    v=Box(lo,hi)
    d2=Box(xlo*xlo,xhi*xhi)
    L=first_log_weight-ONE/(TWO*v)+d2/(TWO*v.pow(2))
    Lp=second_log_weight+ONE/(TWO*v.pow(2))-d2/v.pow(3)
    g=gaussian_box(lo,hi,xlo,xhi)
    mid_value=w_mid*Box(gaussian(mid,xhi).a,gaussian(mid,xlo).b)
    second_derivative=w*g*(squared(L)+Lp)
    return mid_value,unsign(second_derivative)


def certificate(j,k,xlo,xhi,depth):
    assert j>0 and k>0 and D(2)<xlo<=xhi<D(10)
    assert 3<=depth<=20 and all(x.as_tuple().exponent>=-9 for x in (xlo,xhi))
    step=D(1)/(2**depth)
    cutoff=D('0.125')
    crit=Box(j+k)/Box(j)
    peak=(Box(j)/Box(j+k)).pow(j)*(Box(k)/Box(j+k)).pow(k)
    assert Box(cutoff).exp().b<crit.a and xlo*xlo>cutoff
    initial=Box(0,(Box(cutoff)*weight(cutoff,j,k)*gaussian(cutoff,xlo)).b)
    pieces=[initial]
    errors=[]
    positives=ZERO
    negatives=ZERO
    squared_step=Box(step).pow(3)/Box(24)
    for segment,a,b in ((0,cutoff,JUMP1),(1,JUMP1,JUMP2),(2,JUMP2,D('1.5'))):
        total=ZERO
        pp=ZERO
        nn=ZERO
        ep=D(0)
        en=D(0)
        sum_error=D(0)
        wlo=weight(a,j,k)
        e_lo=Box(a).exp()
        for idx in range(int((b-a)/step)):
            lo=a+idx*step
            hi=lo+step
            mid=(lo+hi)/2
            whi=weight(hi,j,k)
            e_hi=Box(hi).exp()
            w=weight_cell(wlo,whi,e_lo,e_hi,crit,peak)
            w_mid=weight(mid,j,k)
            wlo=whi
            exp_t=Box(e_lo.a,e_hi.b)
            exp_minus=exp_t-ONE
            e_lo=e_hi
            first=Box(-j)+Box(k)/exp_minus
            second=(ZERO-Box(k)*exp_t)/exp_minus.pow(2)
            if segment!=2:
                mid_v,d2=gaussian_term(lo,hi,mid,xlo,xhi,w,w_mid,first,second)
            if segment!=2:
                pp=pp+Box(step)*mid_v
                ep=HI.add(ep,HI.multiply(squared_step.b,d2))
            if segment:
                s=SHIFT1 if segment==1 else SHIFT2
                mid_s,d2s=gaussian_term(lo-s,hi-s,mid-s,xlo-2,xhi-2,w,w_mid,first,second)
                term=E*mid_s
                err_term=HI.multiply(squared_step.b,HI.multiply(E.b,d2s))
                if segment==1:
                    nn=nn+Box(step)*term
                    en=HI.add(en,err_term)
                    mid_v=mid_v-term
                    d2=HI.add(d2,HI.multiply(E.b,d2s))
                else:
                    pp=pp+Box(step)*term
                    ep=HI.add(ep,err_term)
                    mid_v=term
                    d2=HI.multiply(E.b,d2s)
            total=total+Box(step)*mid_v
            sum_error=HI.add(sum_error,HI.multiply(squared_step.b,d2))
        positives=positives+pp+Box(ep.copy_negate(),ep)
        negatives=negatives+nn+Box(en.copy_negate(),en)
        pieces.append(total+Box(sum_error.copy_negate(),sum_error))
        errors.append((Box(sum_error)/peak).json())
    result=sum(pieces,ZERO)
    p_bound=(positives/peak).json()
    n_bound=(negatives/peak).json()
    return dict(j=j,k=k,total_order=j+k,spatial_box=[str(xlo),str(xhi)],time_step=str(step),retained_positive_scaled=p_bound,negative_scaled=n_bound,
                normalized_stripped_density=(result/peak).json(),raw_stripped_density=result.json(),
                time_slab_error_bounds=errors,normalized_time_slab_bounds=[(v/peak).json() for v in pieces])


def main():
    p=argparse.ArgumentParser()
    p.add_argument('j',type=int)
    p.add_argument('k',type=int)
    p.add_argument('xlo',type=D)
    p.add_argument('xhi',type=D)
    p.add_argument('depth',type=int)
    p.add_argument('out',type=Path)
    a=p.parse_args()
    result=certificate(a.j,a.k,a.xlo,a.xhi,a.depth)
    a.out.write_text(json.dumps(result,indent=2)+'\n')
    print(result['normalized_stripped_density'])
    assert D(result['normalized_stripped_density'][1])<0

if __name__=='__main__':
    main()
