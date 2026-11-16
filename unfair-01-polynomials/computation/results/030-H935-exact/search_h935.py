"""80 necessary H9 layouts: fixed-point outward interval search, not verifier."""
import json,sys,hashlib
from pathlib import Path
from fractions import Fraction as F
from math import isqrt
S=1<<96
Z=(0,0);ONE=(S,S)
def box(lo,hi):return (lo.numerator*S//lo.denominator,-((-hi.numerator*S)//hi.denominator))
def add(a,b):return(a[0]+b[0],a[1]+b[1])
def neg(a):return(-a[1],-a[0])
def mul(a,b):
    vals=[x*y for x in a for y in b]
    return(min(vals)//S,-((-max(vals))//S))
def pmul(a,b):
    c=[Z]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):c[i+j]=add(c[i+j],mul(x,y))
    return c
def lift(a,p):
    c=[Z]*((len(a)-1)*p+1)
    for i,x in enumerate(a):c[i*p]=x
    return c
def coefficients(data):
    mothers=data['mothers'];blocks=[]
    for i in range(3):
        d=mothers[i]['disk'];r=F(d['radius']);x=F(d['center']['re']);y=F(d['center']['im'])
        re=box(x-r,x+r);im=box(y-r,y+r)
        if i==0:blocks.append([neg(re),ONE])
        else:blocks.append([add(mul(re,re),mul(im,im)),neg(add(re,re)),ONE])
    u=isqrt(5*S*S);sqrt5=(u,u+1)
    # coefficients -t_j=(1 +/- sqrt(5))/2; j=1 uses minus.
    cyc=[[ONE,((S-sqrt5[1])//2,-((-(S-sqrt5[0]))//2)),ONE],
         [ONE,((S+sqrt5[0])//2,-((-(S+sqrt5[1]))//2)),ONE]]
    return blocks,cyc
def main():
    source=Path(sys.argv[1]);data=json.loads(source.read_text());g,c=coefficients(data);rows=[]
    for e in (0,1):
        for choice in (range(1) if e==0 else range(4)):
            for cm in range(16):
                A=[ONE];B=[ONE]
                for i,q in enumerate(g):
                    for p in (3,5):
                        selected=e==1 and p==(3 if i==0 or not(choice>>(i-1)&1) else 5)
                        if selected:A=pmul(A,lift(q,p))
                        else:B=pmul(B,lift(q,p))
                for i,q in enumerate(c):
                    for p in (3,5):
                        selected=cm>>(i+(0 if p==3 else 2))&1
                        if selected:A=pmul(A,lift(q,p))
                        else:B=pmul(B,lift(q,p))
                candidates=[(v[1],side,k) for side,P in enumerate((A,B)) for k,v in enumerate(P) if v[1]<0]
                row={'e':e,'choice':choice,'cyclo':cm}
                if candidates:
                    upper,side,k=min(candidates);row.update(status='negative',side=side,degree=k,upper_numerator=upper,denominator=S)
                else:row['status']='unresolved'
                rows.append(row)
    result={'format':'H935-fixed-point-search-v1','mother_input_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'scale':S,'rows':rows}
    Path(sys.argv[2]).write_text(json.dumps(result,indent=2)+'\n')
    print({'total':len(rows),'negative':sum(r['status']=='negative' for r in rows),'remaining':[r for r in rows if r['status']!='negative']})
if __name__=='__main__':main()
