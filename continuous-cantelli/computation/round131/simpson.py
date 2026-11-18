"""Strict positivity on spatial unions from a fourth-order interval Simpson bound."""
import argparse,json,os,time
from functools import lru_cache
from pathlib import Path
from certify import D,Box,LO,HI,ZERO,ONE,TWO,E,gaussian,gaussian_box,weight
from certify_midpoint import squared,unsign,weight_cell

CUTOFF=D('.375')
DIAGONAL=948
SPACE=D('0.000001')
HALF=D('.5')
END=D('1.5')
THREE,FOUR,SIX,TWELVE=map(Box,(3,4,6,12))

@lru_cache(maxsize=12)
def cells(j,k,depth):
    h=D(1)/(2**depth)
    peak=(Box(j)/Box(j+k)).pow(j)*(Box(k)/Box(j+k)).pow(k)
    crit=Box(j+k)/Box(j)
    factor=(Box(h).pow(5)/Box(2880)).b
    ans=[]
    a=CUTOFF
    wl=weight(a,j,k)
    el=Box(a).exp()
    while a<END:
        b=a+h;m=(a+b)/2
        wm=weight(m,j,k);wr=weight(b,j,k);er=Box(b).exp()
        w=weight_cell(wl,wr,el,er,crit,peak)/peak
        averages=(wl/peak*Box(h)/SIX,wm/peak*Box(h)*FOUR/SIX,wr/peak*Box(h)/SIX)
        et=Box(el.a,er.b);B=et-ONE
        L1=Box(-j)+Box(k)/B
        L2=(ZERO-Box(k)*et)/B.pow(2)
        L3=Box(k)*et*(et+ONE)/B.pow(3)
        L4=(ZERO-Box(k)*et*(et.pow(2)+FOUR*et+ONE))/B.pow(4)
        ranges=[]
        for shift in (D(0),D('.375'),D('.875')):
            if a<=shift:
                ranges.append(None);continue
            va,vm,vb=a-shift,m-shift,b-shift
            v=Box(va,vb)
            inv=ONE/v
            q2=inv.pow(2);q3=q2*inv;q4=q3*inv;q5=q4*inv
            constants=(L1-inv/TWO,L2+q2/TWO,L3-q3,L4+THREE*q4)
            coefficients=(q2/TWO,ZERO-q3,THREE*q4,ZERO-TWELVE*q5)
            ranges.append((va,vm,vb,constants,coefficients))
        terms=(0,) if a<HALF else ((0,1) if a<1 else (2,))
        ans.append((averages,w,[ranges[i] for i in terms],terms,factor))
        a,wl,el=b,wr,er
    return ans


def integrals(j,k,y,depth):
    p=ZERO;n=ZERO
    ep=D(0);en=D(0)
    for weights,w,segments,terms,factor in cells(j,k,depth):
        for data,index in zip(segments,terms):
            va,vm,vb,constants,coefficients=data
            x=y if index==0 else y-2
            square=Box(x*x)
            A,B,C,F=(c+b*square for c,b in zip(constants,coefficients))
            a2=squared(A)
            fourth=a2.pow(2)+SIX*a2*B+THREE*squared(B)+FOUR*A*C+F
            bound=gaussian_box(va,vb,x,x)*w*fourth
            err=HI.multiply(factor,unsign(bound))
            v=sum((ww*gaussian(tt,x) for ww,tt in zip(weights,(va,vm,vb))),ZERO)
            if index!=0:
                v=E*v
                err=HI.multiply(E.b,err)
            if index==1:
                n=n+v;en=HI.add(en,err)
            else:
                p=p+v;ep=HI.add(ep,err)
    return p+Box(ep.copy_negate(),ep),n+Box(en.copy_negate(),en)


def endpoint(j,x):
    y=D(x)*SPACE
    result=D('-Infinity')
    for depth in (8,9,10):
        p,n=integrals(j,DIAGONAL-j,y,depth)
        if p.a>0:
            result=LO.subtract(1,HI.divide(n.b,p.a))
            if result>D('0.000012'):
                return result,depth,p.json(),n.json()
    return result,depth,p.json(),n.json()


def verify(j,intervals,deadline):
    found=[];failed=[];endpoints={}
    stack=[]
    for a,b in intervals:
        lo=int(LO.divide(D(a),SPACE).to_integral_value(rounding='ROUND_FLOOR'))
        hi=int(HI.divide(D(b),SPACE).to_integral_value(rounding='ROUND_CEILING'))
        stack.append((lo,hi))
    while stack:
        a,b=stack.pop()
        if time.monotonic()>deadline:
            stack.append((a,b));failed.extend(stack);break
        for x in (a,b):
            if x not in endpoints:
                endpoints[x]=endpoint(j,x)
        h=D(b-a)*SPACE
        loss=D('5.5')*h*h   # F'' <= 44, interpolate with loss 44h²/8.
        bound=LO.subtract(min(endpoints[a][0],endpoints[b][0]),loss)
        if bound>0:
            found.append([a,b,str(bound)])
        elif b-a<=125:
            failed.append((a,b))
        else:
            m=(a+b)//2
            stack.extend(((m,b),(a,m)))
    return dict(j=j,covered=found,unresolved=failed,endpoints={x:[str(b),d,p,n] for x,(b,d,p,n) in endpoints.items()})


def main():
    p=argparse.ArgumentParser();p.add_argument('shard',type=int);p.add_argument('out',type=Path);p.add_argument('--limit',type=float,default=108);a=p.parse_args()
    rows=json.loads(Path('remaining.json').read_text())
    selected=[(int(j),r) for j,r in rows.items() if int(j)%15==a.shard]
    # Prioritize the weakest cases close to 551, then the rest.
    selected.sort(key=lambda r:abs(r[0]-551))
    start=time.time();cpu=time.process_time();deadline=time.monotonic()+a.limit
    result=[]
    for j,remaining in selected:
        if time.monotonic()>deadline:break
        result.append(verify(j,remaining,deadline))
    a.out.write_text(json.dumps(dict(shard=a.shard,pid=os.getpid(),start=start,end=time.time(),cpu=time.process_time()-cpu,results=result),indent=2)+'\n')

if __name__=='__main__': main()
