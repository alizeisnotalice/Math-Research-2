#!/usr/bin/env python3
"""Exact rectangle-union upper bound for a 2D affine mass family.

Input has base, copy_factors, delta, a_bounds, b_bounds and corners ordered
(a0,b0),(a1,b0),(a0,b1),(a1,b1). Each corner is a complete periodic input.
Labels follow (copy_x,copy_y,base_row) order. Four corners establish affine
mass, nonnegativity, mass conservation, and the uniform dominated tail.
"""
import argparse,bisect,hashlib,json,math,platform,sys,time
from fractions import Fraction as F
from pathlib import Path
from periodic_halo import union_area


class MassFamilyEnvelope:
    def __init__(self, data, bits=40, R=F(5)):
        self.bits,self.R=bits,F(R)
        if type(bits) is not int or bits<0 or self.R<=0:
            raise ValueError('invalid bits or cutoff')
        self.a0,self.a1=map(F,data['a_bounds']);self.b0,self.b1=map(F,data['b_bounds'])
        if not (self.a0<self.a1 and self.b0<self.b1):raise ValueError('nondegenerate family rectangle required')
        factors=data['copy_factors']
        if len(factors)!=2 or any(type(i) is not int or i<=0 for i in factors):raise ValueError('invalid copy factors')
        nxcopy,nycopy=factors;delta=F(data['delta'])
        base=data['base'];self.base=base;self.copy_factors=tuple(factors);self.delta=delta;bx,by,lam=(F(base[k]) for k in ('Px','Py','lam'))
        if min(bx,by,lam)<=0 or delta<0:raise ValueError('invalid base or delta')
        base_pts=[(F(x)%bx,F(y)%by,F(w)) for x,y,w in base['pts']]
        if not base_pts or any(w<=0 for _,_,w in base_pts):raise ValueError('base mass must be positive per label')
        corners=data['corners']
        if len(corners)!=4:raise ValueError('exactly four corners required')
        px,py=bx*nxcopy,by*nycopy;self.px,self.py,self.lam=px,py,lam
        if any((F(c['Px']),F(c['Py']),F(c['lam']))!=(px,py,lam) for c in corners):raise ValueError('period or level mismatch')
        n=len(base_pts);size=n*nxcopy*nycopy
        if any(len(c['pts'])!=size for c in corners):raise ValueError('label coverage mismatch')
        rows=[];maxdist=F(0);kappa=F(0)
        for a in range(nxcopy):
          for b in range(nycopy):
            for i,(x0,y0,w0) in enumerate(base_pts):
              idx=(a*nycopy+b)*n+i
              points=[(F(c['pts'][idx][0])%px,F(c['pts'][idx][1])%py,F(c['pts'][idx][2])) for c in corners]
              x,y=points[0][:2];weights=[p[2] for p in points]
              if any(p[:2]!=(x,y) for p in points) or min(weights)<0:raise ValueError('geometry differs or negative corner mass')
              if weights[0]+weights[3]!=weights[1]+weights[2]:raise ValueError('mass is not affine')
              dx=(x-(x0+a*bx))%px;dy=(y-(y0+b*by))%py
              distance=max(min(dx,px-dx),min(dy,py-dy));maxdist=max(maxdist,distance)
              if distance>delta:raise ValueError('displacement exceeds delta')
              kappa=max(kappa,*(w/w0 for w in weights))
              rows.append((x,y,*weights))
        totals=[sum((row[2+i] for row in rows),F(0)) for i in range(4)]
        if totals[0]<=0 or len(set(totals))!=1:raise ValueError('total mass is not positive and constant')
        self.mass=totals[0];self.points=rows;base_mass=sum((p[2] for p in base_pts),F(0))
        self.tail_bound=kappa*base_mass/(bx*by)*(1+(bx+2*delta)/(2*self.R))*(1+(by+2*delta)/(2*self.R))
        if self.tail_bound>=lam:raise ValueError('uniform dominated tail cutoff fails')
        D=1<<bits;S=math.lcm(D,px.denominator,py.denominator,self.R.denominator,*(v.denominator for row in rows for v in row[:2])); W=math.lcm(*(v.denominator for row in rows for v in row[2:])); self.scale,self.weight_scale=S,W
        X,Y,Ri=int(px*S),int(py*S),int(self.R*S);atoms={}
        for x,y,w00,w10,w01,w11 in rows:
            xi,yi=int(x*S),int(y*S); ws=[int(v*W) for v in (w00,w10,w01,w11)]
            for ix in range((-Ri-xi)//X,(X+Ri-xi)//X+1):
              xx=xi+ix*X
              if not -Ri<=xx<=X+Ri:continue
              for iy in range((-Ri-yi)//Y,(Y+Ri-yi)//Y+1):
                yy=yi+iy*Y
                if -Ri<=yy<=Y+Ri:
                    acc=atoms.setdefault((xx,yy),[0,0,0,0])
                    for k in range(4):acc[k]+=ws[k]
        xs=sorted({x for x,y in atoms});ys=sorted({y for x,y in atoms});nx,ny=len(xs),len(ys)
        grids=[[[0]*(ny+1) for _ in range(nx+1)] for _ in range(3)]
        xm={x:i+1 for i,x in enumerate(xs)};ym={y:i+1 for i,y in enumerate(ys)}
        for (x,y),v in atoms.items():
            for p,val in zip(grids,(v[0],v[1]-v[0],v[2]-v[0])):p[xm[x]][ym[y]]+=val
        for p in grids:
          for i in range(1,nx+1):
            acc=0
            for j in range(1,ny+1):acc+=p[i][j];p[i][j]=p[i-1][j]+acc
        def box(p,i,k,j,l):return p[k+1][l+1]-p[i][l+1]-p[k+1][j]+p[i][j]
        A,B=lam.denominator,4*W*lam.numerator; retained=[]; enumerated=0
        for i,a in enumerate(xs):
          for k in range(i,bisect.bisect_right(xs,a+2*Ri)):
            b=xs[k]
            if b-Ri>=X or a+Ri<=0:continue
            for j,c in enumerate(ys):
              for l in range(j,bisect.bisect_right(ys,c+2*Ri)):
                d=ys[l]
                if d-Ri>=Y or c+Ri<=0:continue
                enumerated+=1;m,da,db=(box(p,i,k,j,l) for p in grids);maxm=max(m,m+da,m+db,m+da+db)
                if maxm<=0:continue
                span=max(b-a,d-c)
                if span*span*B>=4*S*S*maxm*A:continue
                root=math.isqrt(maxm*A*D*D//B);hi=min(Ri,(root+1)*(S//D))
                if max(0,b-hi)>=min(X,a+hi) or max(0,d-hi)>=min(Y,c+hi):continue
                retained.append((a,b,c,d,m,da,db))
        self.boxes=retained
        minimum_corner_mass=min(v for row in rows for v in row[2:])
        self.meta={'bits':bits,'dyadic_outward_radius_error_at_most':str(F(1,D)),'radius':str(self.R),'mass':str(self.mass),'minimum_corner_label_mass':str(minimum_corner_mass),'all_corner_labels_strictly_positive':minimum_corner_mass>0,'uniform_kappa':str(kappa),'delta':str(delta),'actual_maximum_distance':str(maxdist),'uniform_tail_bound':str(self.tail_bound),'tail_margin':str(lam-self.tail_bound),'window_atoms_distinct':len(atoms),'x_coordinates':nx,'y_coordinates':ny,'enumerated_boxes':enumerated,'retained_boxes':len(retained)}

    def point(self,a,b):
        a,b=F(a),F(b)
        if not (self.a0<=a<=self.a1 and self.b0<=b<=self.b1):raise ValueError('point outside family')
        u=(a-self.a0)/(self.a1-self.a0);v=(b-self.b0)/(self.b1-self.b0)
        pts=[[str(x),str(y),str(w00+u*(w10-w00)+v*(w01-w00))] for x,y,w00,w10,w01,w11 in self.points]
        return {'Px':str(self.px),'Py':str(self.py),'lam':str(self.lam),'pts':pts}

    def upper(self,a_lo,a_hi,b_lo,b_hi):
        a_lo,a_hi,b_lo,b_hi=map(F,(a_lo,a_hi,b_lo,b_hi))
        if not (self.a0<=a_lo<=a_hi<=self.a1 and self.b0<=b_lo<=b_hi<=self.b1):raise ValueError('cell outside family')
        aa=[(a_lo-self.a0)/(self.a1-self.a0),(a_hi-self.a0)/(self.a1-self.a0)]
        bb=[(b_lo-self.b0)/(self.b1-self.b0),(b_hi-self.b0)/(self.b1-self.b0)]
        Q=math.lcm(*(v.denominator for v in aa+bb));aa=[int(v*Q) for v in aa];bb=[int(v*Q) for v in bb]
        S,W=self.scale,self.weight_scale;D=1<<self.bits;Ri=int(self.R*S);X,Y=int(self.px*S),int(self.py*S);A,B=self.lam.denominator,4*W*Q*self.lam.numerator
        rects=set()
        for a,b,c,d,m,da,db in self.boxes:
            M=m*Q+da*aa[1 if da>=0 else 0]+db*bb[1 if db>=0 else 0]
            if M<=0:continue
            span=max(b-a,d-c)
            if span*span*B>=4*S*S*M*A:continue
            n=M*A*D*D;k=math.isqrt(n//B)
            assert k*k*B<=n<(k+1)*(k+1)*B
            s=min(Ri,(k+1)*(S//D));x0,x1=max(0,b-s),min(X,a+s);y0,y1=max(0,d-s),min(Y,c+s)
            if x0<x1 and y0<y1:rects.add((x0,x1,y0,y1))
        area=F(union_area(rects),S*S)
        return {'upper_area':str(area),'upper_ratio':str(self.lam*area/self.mass),'rectangles_unique':len(rects),'cell':[str(a_lo),str(a_hi),str(b_lo),str(b_hi)]}


def main():
    p=argparse.ArgumentParser();p.add_argument('input');p.add_argument('--cell',nargs=4);p.add_argument('--bits',type=int,default=40);p.add_argument('--radius',default='5');p.add_argument('--output',required=True);args=p.parse_args()
    start=time.time();path=Path(args.input);f=MassFamilyEnvelope(json.loads(path.read_text()),args.bits,F(args.radius));cell=args.cell or [str(f.a0),str(f.a1),str(f.b0),str(f.b1)]
    result={'input_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'python':sys.version,'platform':platform.platform(),'family':f.meta,'bound':f.upper(*cell),'elapsed_seconds':time.time()-start}
    Path(args.output).write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
