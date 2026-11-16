"""Independent outward fixed-point search of the 17 complete-fiber classes."""
import json,sys,hashlib
from fractions import Fraction as F
from pathlib import Path
S=1<<96;Z=(0,0);O=(S,S)
def box(lo,hi):return(lo.numerator*S//lo.denominator,-((-hi.numerator*S)//hi.denominator))
def add(a,b):return(a[0]+b[0],a[1]+b[1])
def neg(a):return(-a[1],-a[0])
def mul(a,b):
    p=[x*y for x in a for y in b];return(min(p)//S,-((-max(p))//S))
def pmul(a,b):
    c=[Z]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):c[i+j]=add(c[i+j],mul(x,y))
    return c
def lift(a,p):
    c=[Z]*((len(a)-1)*p+1)
    for i,x in enumerate(a):c[i*p]=x
    return c
def main():
    path=Path(sys.argv[1]);data=json.loads(path.read_text());blocks=[]
    for mother in sorted(data['mothers'],key=lambda x:x['block']):
        d=mother['disk'];r=F(d['radius']);x=F(d['center']['re']);y=F(d['center']['im'])
        re=box(x-r,x+r);im=box(y-r,y+r)
        blocks.append([add(mul(re,re),mul(im,im)),neg(add(re,re)),O])
    assert len(blocks)==5
    rows=[{'kind':'constant','status':'exact_constant'}]
    for mask in range(16):
        A=[O];B=[O]
        for i,block in enumerate(blocks):
            p=3 if i==0 or not(mask>>(i-1)&1) else 5
            A=pmul(A,lift(block,p));B=pmul(B,lift(block,8-p))
        candidates=[(v[1],side,k) for side,P in enumerate((A,B)) for k,v in enumerate(P) if v[1]<0]
        row={'kind':'mixed_fibers','mask':mask}
        if candidates:
            u,side,k=min(candidates);row.update(status='negative',side=side,degree=k,upper_numerator=u,denominator=S)
        else:row['status']='unresolved'
        rows.append(row)
    out={'format':'H1035-fixed-point-search-v1','mother_source_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'rows':rows}
    Path(sys.argv[2]).write_text(json.dumps(out,indent=2)+'\n');print({'classes':len(rows),'negative':sum(x['status']=='negative' for x in rows),'remaining':[x for x in rows if x['status']!='negative']})
if __name__=='__main__':main()
