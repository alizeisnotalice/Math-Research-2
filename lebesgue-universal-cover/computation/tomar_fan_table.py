"""Build integer-coordinate witnesses at every possible central angle of the frozen 11D tree."""
import argparse,json,hashlib
from fractions import Fraction
from pathlib import Path
import numpy as np
from flint import arb,ctx
p=argparse.ArgumentParser(description=__doc__)
p.add_argument("certdir",type=Path);p.add_argument("out",type=Path);p.add_argument("--subdiv",type=int,default=8)
a=p.parse_args()
h=json.loads((a.certdir/"part0000.json").read_bytes());ctx.prec=128
q=2**28;grid=1024
counts=[96]+[n*a.subdiv for n in (3,5,7,9)]
out=np.zeros((grid+1,sum(counts),2),dtype=np.int32);offset=0
def from_frac(x):
    f=Fraction(x);return arb(f.numerator)/f.denominator
for body,(n,count) in enumerate(zip([0,3,5,7,9],counts)):
    pi=arb.pi();r=arb(1)/2 if body==0 else 1/(2*(pi/(2*n)).cos())
    def vec(phi): return (phi.cos(),phi.sin())
    base=[];vertices=[]
    if body: vertices=[tuple(r*x for x in vec(2*pi*i/n)) for i in range(n)]
    for j in range(count):
        d=vec(2*pi*j/count)
        if body:
            k=int((Fraction(j*n,count)+Fraction(n,2)+Fraction(1,2)).__floor__())%n
            dot=d[0]*vertices[k][0]+d[1]*vertices[k][1]
            rho=dot+(1-r*r+dot*dot).sqrt()
        else: rho=r
        base.append(tuple(value*rho*(1-arb(1)/10**6) for value in d))
    hi=0 if body<2 else h["hi"][3*body-4]
    for index in ([0] if body<2 else range(grid+1)):
        phi=from_frac(hi)*index/grid; cs,sn=vec(phi)
        def rot(v):return (cs*v[0]-sn*v[1],sn*v[0]+cs*v[1])
        rotated=[rot(v) for v in vertices]
        for j,xy in enumerate(base):
            val=rot(xy);nx,ny=(round(float(x.mid())*q) for x in val)
            x,y=arb(nx)/q,arb(ny)/q
            if body:
                assert all((x-px)**2+(y-py)**2<1 for px,py in rotated)
            else: assert x*x+y*y<arb(1)/4
            out[index,offset+j]=[nx,ny]
    if body<2: out[1:,offset:offset+count]=out[0,offset:offset+count]
    offset+=count
assert (np.abs(out.astype(np.int64))<q).all()
fine=2**44
orig_lo=[Fraction(x)*fine for x in h["lo"]];orig_hi=[Fraction(x)*fine for x in h["hi"]]
assert all(x.denominator==1 for x in orig_lo+orig_hi)
np.savez_compressed(a.out,points=out,counts=np.array(counts,dtype=np.int32),scale=np.array(q),grid=np.array(grid),
    original_lo=np.array(list(map(int,orig_lo)),dtype=np.int64),original_hi=np.array(list(map(int,orig_hi)),dtype=np.int64),
    original_scale=np.array(fine),angle_indices=np.array([2,5,8],dtype=np.int32),
    manifest_sha256=np.array(hashlib.sha256((a.certdir/"SHA256SUMS").read_bytes()).hexdigest()),
    generator_sha256=np.array(hashlib.sha256(Path(__file__).read_bytes()).hexdigest()))
print(json.dumps({"table":str(a.out),"shapes":counts,"scale":q,"angular_grid":grid,"all_witnesses_strictly_inside":True}))
