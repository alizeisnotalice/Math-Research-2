"""Exact all-orientation non-fitting into an area illumination body around an 8D box."""
from argparse import ArgumentParser
from fractions import Fraction as Q
from hashlib import sha256
import json
from math import floor,gcd

from flint import arb,ctx
import numpy as np
from scipy.optimize import linprog

from certify_rotation_inballs import hull

# Fixed decimal centre from Tomar v9 independent/opt8_point.json, rounded to 10 decimals.
CENTER=tuple(map(Q,('-0.0025482897','-0.0019350278','0.0287183302','0.0250869718','-0.0032906511','0.7016242528','0.0028934171','-0.0127181545')))
GRID=10**10
TAU=Q(41831,50000) # 0.83662
HALF=Q(1,50000)

def ball(q):
    q=Q(q);return arb(q.numerator)/q.denominator

def det(a,b):return a[0]*b[1]-a[1]*b[0]

def inner_polygon(half=HALF):
    p=list(map(ball,CENTER));pi=arb.pi();points=[]
    def emit(center,radius,a,step,count):
        for j in range(count):
            sn,cs=(a+j*step).sin_cos();pt=(center[0]+radius*cs,center[1]+radius*sn);out=[]
            for q in pt:
                n=round(float(q.mid())*GRID);assert abs(q*GRID-n)<1;out.append(n)
            points.append(tuple(out))
    K=64;emit((0,0),arb(1)/2,arb(0),pi/(4*K),8*K)
    for n,rot,x,y in [(3,arb(0),p[0],p[1]),(5,p[2],p[3],p[4]),(7,p[5],p[6],p[7])]:
        rn=1/(2*(pi/(2*n)).cos());assert rn<arb(3)/5
        for j in range(n):
            ph=rot+2*pi*j/n;sn,cs=ph.sin_cos();emit((x+rn*cs,y+rn*sn),arb(1),ph+pi-pi/(2*n),pi/(n*K),K+1)
    # Rotated regular Reuleaux bodies lie in their radius-Rn circumdisks (also verified analytically in proof).
    eta=Q(2,GRID);delta=(Q(99,70)+Q(3,5))*half;alpha=(Q(1,2)-delta)/(Q(1,2)+eta)
    assert Q(99,70)**2>2
    return hull(points),alpha

def cuts_for(p,alpha):
    # P=alpha/grid*p. Build necessary halfplanes from arbitrary subsets of edges.
    s=alpha/GRID;keyed={};D=np.array(p,dtype=float)*float(s);edges=np.roll(D,-1,axis=0)-D;norms=np.c_[edges[:,1],-edges[:,0]];bf=(norms*D).sum(axis=1)
    area=sum(det(a,b) for a,b in zip(p,p[1:]+p[:1]))*s*s/2;budget=TAU-area
    assert budget>0
    for phi in np.arange(256)*np.pi/128:
        direction=np.array([np.cos(phi),np.sin(phi)]);lo,hi=.45,.65
        for _ in range(38):
            mid=(lo+hi)/2
            if np.maximum(0,norms@(mid*direction)-bf).sum()>2*float(budget):hi=mid
            else:lo=mid
        inds=tuple(map(int,np.where(norms@(hi*direction)>bf)[0]))
        if not inds or inds in keyed:continue
        nx=ny=intercept=0
        for j in inds:
            a,b=p[j],p[(j+1)%len(p)];nx+=b[1]-a[1];ny+=a[0]-b[0];intercept+=det(a,b)
        div=gcd(abs(nx),abs(ny));assert div>0
        keyed[inds]=((nx//div,ny//div),(s*intercept+2*budget/s)/div)
    return list(keyed.values()),area

def material(n):
    # Finite lattice points inside B9 by checking all nine defining unit disks.
    pi=arb.pi();R=1/(2*(pi/(2*n)).cos());centers=[];out=[]
    for k in range(n):
        s,c=(2*pi*k/n).sin_cos();centers.append((R*c,R*s))
    for k,(x,y) in enumerate(centers):
        for j in range(17):
            ph=2*pi*k/n+pi-pi/(2*n)+j*pi/(16*n);s,c=ph.sin_cos();v=(x+c,y+s)
            q=tuple(round(float(z.mid())*GRID*(1-1e-7)) for z in v)
            a,b=[arb(z)/GRID for z in q]
            assert all((a-u)**2+(b-vv)**2<1 for u,vv in centers)
            out.append(q)
    return hull(out)

def run(M):
    ctx.prec=112;pi=arb.pi();p,alpha=inner_polygon();cuts,area=cuts_for(p,alpha)
    rows,least,Q9=certify_rows(cuts,M)
    return {'tau':str(TAU),'center':[str(q) for q in CENTER],'box_halfwidth':str(HALF),'alpha':str(alpha),'inner_vertices':len(p),'inner_area':str(area),
            'cuts':len(cuts),'integer_material_B9':len(Q9),'angular_cells':M,'full_angle':'[0,2*pi/9]','minimum_dual_gap_lower':str(Q(least,10**12)),
            'rows':rows,'scope':'strictly excludes the indicated 8D base box at tau for every B9 translation and rotation; global remaining boxes unresolved'}

def certify_rows(cuts,M):
    pi=arb.pi();Q9=material(9);fm=np.array([q for q,b in cuts],dtype=float);ln=np.linalg.norm(fm,axis=1);F=fm/ln[:,None];fb=np.array([float(b) for q,b in cuts])/ln
    vertices=np.array(Q9,dtype=float)/GRID;rows=[];least=None
    for j in range(M):
        angle=(j+.5)*2*np.pi/(9*M);cs,sn=np.cos(angle),np.sin(angle)
        verts=vertices@np.array([[cs,sn],[-sn,cs]]);mat=F@verts.T;ids=mat.argmax(axis=1);sup=mat.max(axis=1)
        res=linprog([0,0,1],A_ub=np.c_[F,-np.ones(len(F))],b_ub=fb-sup,bounds=[(None,None)]*3,method='highs');assert res.success
        active=np.argsort(abs(res.ineqlin.marginals))[-3:].tolist();normals=[cuts[i][0] for i in active]
        ws=[det(normals[1],normals[2]),det(normals[2],normals[0]),det(normals[0],normals[1])]
        if sum(ws)<0:ws=[-w for w in ws]
        assert min(ws)>=0 and sum(ws)>0 and all(sum(w*q[k] for w,q in zip(ws,normals))==0 for k in range(2))
        radius=pi/(9*M);cang=(2*j+1)*radius;sn,cs=cang.sin_cos();gap=arb(0)
        for w,i in zip(ws,active):
            q=np.asarray(Q9[ids[i]],dtype=object);x,y=[arb(int(z))/GRID for z in q];nx,ny=cuts[i][0]
            vx=cs*x-sn*y;vy=sn*x+cs*y;dot=nx*vx+ny*vy;der=-nx*vy+ny*vx
            rem=(arb(nx*nx+ny*ny)*(x*x+y*y)).sqrt()*radius**2/2
            gap+=w*(dot-abs(der)*radius-rem-ball(cuts[i][1]))
        # Normalize solely for small readable scalar comparisons; exact balance is retained above.
        denom=sum(w*(abs(q[0])+abs(q[1])) for w,q in zip(ws,normals))
        margin=gap/denom;lower=floor(float(margin.lower())*10**12)-1
        assert margin>arb(lower)/10**12>0,(j,float(margin.mid()))
        least=lower if least is None else min(least,lower)
        rows.append({'active':active,'vertices':[int(ids[i]) for i in active],'dual_lower_1e12':lower})
    return rows,least,Q9

if __name__=='__main__':
    ap=ArgumentParser();ap.add_argument('--cells',type=int,default=128);a=ap.parse_args();print(json.dumps(run(a.cells),separators=(',',':')))
