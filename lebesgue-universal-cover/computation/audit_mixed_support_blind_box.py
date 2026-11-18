"""Upper-bound every common-inner-polygon mixed-support budget on one preselected full original box."""
import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
import sys
import numpy as np
from flint import arb,ctx

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
from tomar_leaf_box import leaf_box
S=2**44
D=10**8

def sh(p):return sha256(p.read_bytes()).hexdigest()
def q(f):
    f=F(f);return arb(f.numerator)/f.denominator
def det(a,b):return a[0]*b[1]-a[1]*b[0]
def dot(a,b):return a[0]*b[0]+a[1]*b[1]
def ceil_rat(z):
    return F(int((z.upper()*D).ceil().unique_fmpz()),D)

def vertices(n,angle):
    r=1/(2*(arb.pi()/(2*n)).cos())
    out=[]
    for k in range(n):
        si,co=(angle+2*arb.pi()*k/n).sin_cos()
        out.append((r*co,r*si))
    return out

def support_upper(v,n,l):
    vals=[dot(z,n) for z in v]
    mh=max(z.upper() for z in vals)
    ml=max(z.lower() for z in vals)
    ah=l.upper()+min(z.upper() for z in vals)
    if ml>(l/2).upper():return mh
    if mh<(l/2).lower():return ah
    return max(mh,ah)

def main():
    if not __debug__:raise SystemExit("Assertions required.")
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('main','cert','input'):p.add_argument(name,type=Path)
    p.add_argument('--choose-facets',type=Path)
    a=p.parse_args();ctx.prec=192
    corpus=a.main/'finite-point-lift-083741.json'
    records=json.loads(corpus.read_bytes())['rows']
    pool=sorted((sha256(f'guide-3.3-holdout:{r[0]}:{r[1]}'.encode()).hexdigest(),r)
                for r in records if r[5]>=300000)
    assert len(pool)==35 and pool[0][1][:2]==[1202,32259]
    part,idx=pool[0][1][:2]
    raw=(a.cert/'SHA256SUMS').read_bytes()
    assert sha256(raw).hexdigest()=='56cbe5dfdadc6b5aa07cf6ec9acf4f637a3c4fa78db54df0ae2d78916d292410'
    pins={n.removeprefix('./'):v for v,n in (s.split()for s in raw.decode().splitlines())}
    stem=f'part{part:04}'
    for ext in ('.json','.npy'):
        assert sh(a.cert/(stem+ext))==pins[stem+ext]
    meta=json.loads((a.cert/(stem+'.json')).read_bytes())
    codes=np.load(a.cert/(stem+'.npy'),allow_pickle=False)
    prefix=np.r_[1,1+np.cumsum(np.where(codes<11,1,-1))]
    assert meta['part']==part and codes[idx]>=253 and prefix[-1]==0 and prefix[:-1].min()==1
    c,w=leaf_box(prefix,codes,part,idx,meta)
    lo=[F(x)-F(y)for x,y in zip(c,w)];hi=[F(x)+F(y)for x,y in zip(c,w)]
    assert all(0<=lo[i]<=hi[i] for i in (2,5,8)) and q(hi[2])<arb.pi()/15
    assert all(max(abs(lo[j]),abs(hi[j]))**2+max(abs(lo[j+1]),abs(hi[j+1]))**2<F('.04')**2 for j in (0,3,6,9))
    rec=json.loads(a.input.read_bytes());ns=[tuple(map(int,n))for n in rec['normals']]
    ref=[F(v,S)for v in rec['reference']]
    assert all(l<=v<=h for l,v,h in zip(lo,ref,hi))
    n30=a.main/'finite-point-leaf-input.json'
    assert sh(n30)=='b99b3bf66adf904b8b99ba1e1c7f8e75abc1ea6d8e00728b9553c20b961ffebc'
    data=json.loads(n30.read_bytes());pts=[tuple(F(v,data['denominator'])for v in x)for x in data['points']]
    assert len(pts)==30 and set(pts)=={(x,-y)for x,y in pts}
    angle,tx,ty=map(lambda t:F(t,rec['point_pose_scale']),rec['point_pose'])
    si,co=q(angle).sin_cos()
    pm=[(q(tx)+co*q(x)-si*q(y),q(ty)+si*q(x)+co*q(y))for x,y in pts]
    trial={};reference={}
    for i,n in enumerate((3,5,7,9)):
        k=3*i;opts=[F(0)] if i==0 else [lo[k-1],(lo[k-1]+hi[k-1])/2,hi[k-1]]
        trial[n]=[vertices(n,q(t)) for t in opts]
        reference[n]=vertices(n,q(ref[k-1]) if i else q(0))
    offsets=[];href=[]
    for normal in ns:
        L=q(dot(normal,normal)).sqrt()
        common=[L.upper()/2];host=[L.upper()/2,max(dot(z,normal).upper()for z in pm)]
        for i,n in enumerate((3,5,7,9)):
            k=3*i
            translation=min(normal[0]*lo[k],normal[0]*hi[k])+min(normal[1]*lo[k+1],normal[1]*hi[k+1])
            common.append((q(translation)+min(support_upper(v,normal,L)for v in trial[n])).upper())
            host.append((q(dot(ref[k:k+2],normal))+support_upper(reference[n],normal,L)).upper())
        offsets.append(ceil_rat(max(common)));href.append(ceil_rat(max(host)))
    if a.choose_facets:
        from scipy.spatial import ConvexHull
        dual=np.asarray(ns,dtype=float)/np.asarray(offsets,dtype=float)[:,None]
        ids=ConvexHull(dual).vertices.tolist()
        minimal={k:rec[k] for k in ('reference','point_pose','point_pose_scale')}
        minimal['normals']=[ns[i]for i in ids]
        minimal['holdout']=[part,idx]
        a.choose_facets.write_text(json.dumps(minimal,separators=(',',':'))+'\n')
        return
    # Floating-point geometry is not part of the certificate. Positively ordered
    # rational normals and all rational line-line intersections define Q exactly.
    assert len(ns)==len(set(ns)) and all(z>0 for z in offsets)
    points=[]
    for j,a0 in enumerate(ns):
        b0=ns[(j+1)%len(ns)];z=det(a0,b0);assert z>0
        r=offsets[j];s=offsets[(j+1)%len(ns)]
        vertex=((r*b0[1]-a0[1]*s)/z,(a0[0]*s-r*b0[0])/z)
        assert all(dot(vertex,n)<=v for n,v in zip(ns,offsets))
        points.append(vertex)
    mixed=F(0);area=F(0)
    for i,n in enumerate(ns):
        old,new=points[i-1],points[i];edge=tuple(v-u for u,v in zip(old,new))
        alpha=det(n,edge)/dot(n,n)
        assert alpha>0 and edge==(-alpha*n[1],alpha*n[0])
        mixed+=alpha*href[i]/2
        area+=det(old,new)/2
    assert 0<area<mixed<F('.83741')
    upper=F((mixed*10**12).__ceil__(),10**12)
    print(json.dumps(dict(original=[part,idx],R169_nodes=pool[0][1][3],R169_p30_terminals=pool[0][1][5],
        blind_pool_size=len(pool),blind_hash=pool[0][0],normals=len(ns),point_hand='both by rational x-axis symmetry',
        common_outer_area_upper_at_1e12=str(F((area*10**12).__ceil__(),10**12)),mixed_area_upper_at_1e12=str(upper),
        margin_at_083741=str(F('.83741')-upper),margin_at_083742=str(F('.83742')-upper),
        exact_mixed_upper_sha256=sha256(str(mixed).encode()).hexdigest(),input_sha256=sh(a.input),source_sha256=sh(Path(__file__)),
        R169_ledger_sha256=sh(corpus),pointset_sha256=sh(n30),domain='all placements in the whole original 11D box',
        scope='any single common inner polygon and any balanced shared-capacity rows for the five bodies and the original P30'),indent=2))
if __name__=='__main__':main()
