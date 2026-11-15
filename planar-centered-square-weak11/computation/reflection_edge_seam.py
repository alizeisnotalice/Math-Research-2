#!/usr/bin/env python3
"""Certify the stabilized rightmost seam using a bootstrap-reduced finite catalog."""
import argparse,bisect,hashlib,json,math,os,platform,sys,time
from fractions import Fraction as F
from pathlib import Path

def ceilq(x):return -(-x.numerator//x.denominator)
def upward_root(x,den):
    n=math.isqrt(x.numerator*den*den//x.denominator)
    assert F(n,den)**2<=x<F(n+1,den)**2
    return F(n+1,den)
def calculate(data,bits=64):
    px,py,lam=(F(data[k]) for k in ('Px','Py','lam'))
    pts=[tuple(map(F,p)) for p in data['pts'] if F(p[2])>0]
    assert all(0<=x<px and 0<=y<py for x,y,m in pts)
    mass=sum(m for x,y,m in pts);rho=mass/(px*py)
    old_bound=rho*(1+px/20)*(1+py/20);reflected_bound=rho*(1+px/15)*(1+py/30)
    assert old_bound<lam and reflected_bound<lam
    xmax=max(x for x,y,m in pts);R=F(10);boot=[]
    for step in range(30):
        nc=ceilq(R/px);nr=ceilq(2*R/py)+1
        row=[(i,j,x+i*px,y,m) for i in range(1-nc,1) for j,(x,y,m) in enumerate(pts) if x+i*px>xmax-R]
        rowmass=sum(t[4] for t in row);B=nr*rowmass
        new=upward_root(B/(4*lam),1000)
        boot.append({'R':str(R),'columns':nc,'rows':nr,'row_mass':str(rowmass),'mass_upper':str(B),'next_strict_radius':str(new),'row_labels':[[i,j] for i,j,x,y,m in row]})
        if new>=R:break
        R=new
    else:raise AssertionError('bootstrap failed to stop')
    raw=[(x,y+j*py,m) for i,label,x,y,m in row for j in range(nr)]
    D=1<<bits;S=math.lcm(D,*(z.denominator for p in raw for z in p[:2]));W=math.lcm(*(p[2].denominator for p in raw));Rint=int(R*S)
    atoms={}
    for x,y,m in raw:
        key=(int(x*S),int(y*S));atoms[key]=atoms.get(key,0)+int(m*W)
    xs=sorted({x for x,y in atoms});ys=sorted({y for x,y in atoms});nx,ny=len(xs),len(ys)
    xm={x:i for i,x in enumerate(xs)};ym={y:j for j,y in enumerate(ys)}
    p=[[0]*(ny+1) for _ in range(nx+1)]
    for (x,y),m in atoms.items():p[xm[x]+1][ym[y]+1]+=m
    for i in range(1,nx+1):
        a=0
        for j in range(1,ny+1):a+=p[i][j];p[i][j]=p[i-1][j]+a
    A=lam.denominator;B=4*W*lam.numerator
    best=None;owner=None;count={'boxes':0,'span_rejected':0,'cannot_beat_lower':0,'root_checked':0}
    for (x,y),m in atoms.items():
        rr=math.isqrt(m*A*D*D//B);value=x+rr*(S//D)
        if best is None or value>best:best=value;owner={'bbox':[str(F(x,S)),str(F(x,S)),str(F(y,S)),str(F(y,S))],'mass':str(F(m,W))}
    for i,a in enumerate(xs):
      for k in range(i,bisect.bisect_left(xs,a+2*Rint)):
        b=xs[k]
        for j,c in enumerate(ys):
          for l in range(j,bisect.bisect_left(ys,c+2*Rint)):
            d=ys[l];count['boxes']+=1
            M=p[k+1][l+1]-p[i][l+1]-p[k+1][j]+p[i][j]
            span=max(b-a,d-c)
            if not M or span*span*B>=4*S*S*M*A:count['span_rejected']+=1;continue
            assert best>a
            if S*S*M*A<=B*(best-a)**2:count['cannot_beat_lower']+=1;continue
            rr=math.isqrt(M*A*D*D//B);count['root_checked']+=1
            value=a+rr*(S//D)
            if value>best:best=value;owner={'bbox':[str(F(a,S)),str(F(b,S)),str(F(c,S)),str(F(d,S))],'mass':str(F(M,W))}
    lo=F(best,S);hi=lo+F(1,D)
    assert lo>xmax
    return {'status':'PASS','Px':str(px),'Py':str(py),'lam':str(lam),'mass':str(mass),'density':str(rho),'old_tail_R10':str(old_bound),'reflection_tail_R15':str(reflected_bound),'bootstrap':boot,'final_owner_radius':str(R),'catalog_columns':nc,'catalog_rows':nr,'catalog_atoms':len(atoms),'catalog_x':nx,'catalog_y':ny,'counts':count,'c_edge_lower':str(lo),'c_edge_upper':str(hi),'rational_reflection_axis':str(hi),'c_edge_interval_width':str(hi-lo),'max_lower_witness':owner,'bits':bits,'N_stabilization_sufficient':max(nc,nr),'python':sys.version,'platform':platform.platform(),'scope':'Exact period-count tail, strict owner-radius bootstrap, complete source-bbox maximum enclosure. Upper endpoint is an authorized rational outward reflection axis, not a claim that the true seam is rational.'}
def main():
    ap=argparse.ArgumentParser();ap.add_argument('input');ap.add_argument('--output',required=True);ap.add_argument('--bits',type=int,default=64);args=ap.parse_args()
    assert os.path.ismount('/Volumes/data')
    t=time.time();path=Path(args.input);result=calculate(json.loads(path.read_text()),args.bits);result.update(input_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),elapsed_seconds=time.time()-t)
    out=Path(args.output);out.parent.mkdir(parents=True,exist_ok=True);out.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n');print(json.dumps({k:result[k] for k in ('status','catalog_atoms','catalog_x','catalog_y','counts','c_edge_lower','c_edge_upper','elapsed_seconds')},ensure_ascii=False))
if __name__=='__main__':main()
