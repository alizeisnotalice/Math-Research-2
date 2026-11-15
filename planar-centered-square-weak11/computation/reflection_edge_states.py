#!/usr/bin/env python3
"""Complete local halo enclosures for periodic half-planes and reflected seams."""
import argparse,bisect,hashlib,json,math,os,platform,sys,time
from fractions import Fraction as F
from pathlib import Path
from periodic_halo import union_area

def ceilq(x):return -(-x.numerator//x.denominator)
def source_atoms(pts,px,py,win,R,state,C):
    x0,x1,y0,y1=win;atoms={}
    for reflected in ([False,True] if state=='seam-new' else [False]):
      for a,b,m in pts:
        if reflected:lx,ux=2*C-(x1+R),2*C-(x0-R)
        else:lx,ux=x0-R,x1+R
        imin=ceilq((lx-a)/px);imax=(ux-a)//px
        jmin=ceilq((y0-R-b)/py);jmax=(y1+R-b)//py
        if state=='left':imin=max(imin,0)
        if state=='right':imax=min(imax,-1)
        if state=='bottom':jmin=max(jmin,0)
        if state=='top':jmax=min(jmax,-1)
        if state in ('seam-old','seam-new'):imax=min(imax,0)
        for i in range(imin,imax+1):
          x=a+i*px
          if reflected:x=2*C-x
          for j in range(jmin,jmax+1):
            y=b+j*py;atoms[x,y]=atoms.get((x,y),F(0))+m
    return atoms

def calculate(data,state,seam=None,bits=40):
    px,py,lam=(F(data[k]) for k in ('Px','Py','lam'))
    pts=[tuple(map(F,p)) for p in data['pts'] if F(p[2])>0]
    assert all(0<=x<px and 0<=y<py for x,y,m in pts)
    mass=sum(m for x,y,m in pts);rho=mass/(px*py);C=F(seam['rational_reflection_axis']) if seam else None
    if state in ('left','right'):win=(-2*px,2*px,F(0),py)
    elif state in ('bottom','top'):win=(F(0),px,-2*py,2*py)
    elif state in ('seam-old','seam-new'):
        assert seam and seam['input_sha256'] in (INPUT_HASH,)
        assert C>=F(seam['c_edge_upper'])>max(x for x,y,m in pts)
        win=(C-15,C,F(0),py)
    elif state=='bulk':win=(F(0),px,F(0),py)
    else:raise ValueError(state)
    R=F(15 if state=='seam-new' else 10)
    tail_bound=rho*(1+px/R)*(1+py/(2*R)) if state=='seam-new' else rho*(1+px/(2*R))*(1+py/(2*R))
    assert tail_bound<lam
    atoms=source_atoms(pts,px,py,win,R,state,C)
    D=1<<bits;scale=math.lcm(D,R.denominator,*(v.denominator for v in win),*(z.denominator for p in atoms for z in p));W=math.lcm(*(m.denominator for m in atoms.values()))
    raw={(int(x*scale),int(y*scale)):int(m*W) for (x,y),m in atoms.items()}
    xs=sorted({x for x,y in raw});ys=sorted({y for x,y in raw});nx,ny=len(xs),len(ys)
    xm={x:i+1 for i,x in enumerate(xs)};ym={y:j+1 for j,y in enumerate(ys)}
    p=[[0]*(ny+1) for _ in range(nx+1)]
    for (x,y),m in raw.items():p[xm[x]][ym[y]]+=m
    for i in range(1,nx+1):
        acc=0
        for j in range(1,ny+1):acc+=p[i][j];p[i][j]=p[i-1][j]+acc
    A=lam.denominator;B=4*W*lam.numerator;Ri=int(R*scale);X0,X1,Y0,Y1=(int(v*scale) for v in win)
    inner=set();outer=set();count={'boxes':0,'zero_mass':0,'span_rejected':0,'root_checked':0,'inner_empty':0,'outer_empty':0}
    scale2=scale*scale;D2=D*D;factor=scale//D
    xlimits=[bisect.bisect_right(xs,a+2*Ri) for a in xs];ylimits=[bisect.bisect_right(ys,c+2*Ri) for c in ys]
    for i,a in enumerate(xs):
      for k in range(i,xlimits[i]):
        b=xs[k]
        if b-Ri>=X1 or a+Ri<=X0:continue
        rowhi=p[k+1];rowlo=p[i]
        for j,c in enumerate(ys):
          for l in range(j,ylimits[j]):
            d=ys[l]
            if d-Ri>=Y1 or c+Ri<=Y0:continue
            count['boxes']+=1
            M=rowhi[l+1]-rowlo[l+1]-rowhi[j]+rowlo[j]
            if not M:count['zero_mass']+=1;continue
            span=max(b-a,d-c)
            if span*span*B>=4*scale2*M*A:count['span_rejected']+=1;continue
            rr=math.isqrt(M*A*D2//B);count['root_checked']+=1
            assert rr*rr*B<=M*A*D2<(rr+1)*(rr+1)*B
            lo=min(Ri,max(0,rr-1)*factor);hi=min(Ri,(rr+1)*factor)
            for radius,target,label in ((lo,inner,'inner_empty'),(hi,outer,'outer_empty')):
                lft,rgt=max(X0,b-radius),min(X1,a+radius);bot,top=max(Y0,d-radius),min(Y1,c+radius)
                if lft<rgt and bot<top:target.add((lft,rgt,bot,top))
                else:count[label]+=1
    ai,ao=union_area(inner),union_area(outer)
    return {'status':'PASS','state':state,'Px':str(px),'Py':str(py),'lam':str(lam),'cell_mass':str(mass),'density':str(rho),'window':list(map(str,win)),'tail_radius':str(R),'tail_bound':str(tail_bound),'rational_axis':str(C) if C is not None else None,'bits':bits,'scale':str(scale),'weight_scale':str(W),'window_atoms':len(atoms),'window_x':nx,'window_y':ny,'counts':count,'inner_rectangles_unique':len(inner),'outer_rectangles_unique':len(outer),'inner_area':str(F(ai,scale2)),'outer_area':str(F(ao,scale2)),'inner_area_float':float(F(ai,scale2)),'outer_area_float':float(F(ao,scale2)),'python':sys.version,'platform':platform.platform(),'scope':'Complete source-bbox enumeration for the stated infinite half-plane or reflected-half-plane source restricted to the local receiver window; all radii excluded above the independently certified tail.'}
def main():
    global INPUT_HASH
    ap=argparse.ArgumentParser();ap.add_argument('input');ap.add_argument('--state',required=True,choices=['bulk','left','right','bottom','top','seam-old','seam-new']);ap.add_argument('--seam');ap.add_argument('--bits',type=int,default=40);ap.add_argument('--output',required=True);args=ap.parse_args()
    assert os.path.ismount('/Volumes/data');t=time.time();src=Path(args.input);INPUT_HASH=hashlib.sha256(src.read_bytes()).hexdigest();seam=json.loads(Path(args.seam).read_text()) if args.seam else None
    result=calculate(json.loads(src.read_text()),args.state,seam,args.bits);result.update(input_sha256=INPUT_HASH,script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),union_module_sha256=hashlib.sha256(Path(__file__).with_name('periodic_halo.py').read_bytes()).hexdigest(),seam_sha256=hashlib.sha256(Path(args.seam).read_bytes()).hexdigest() if args.seam else None,elapsed_seconds=time.time()-t)
    out=Path(args.output);out.parent.mkdir(parents=True,exist_ok=True);out.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n');print(json.dumps({k:result[k] for k in ('status','state','window_atoms','window_x','window_y','counts','inner_area','outer_area','elapsed_seconds')},ensure_ascii=False))
if __name__=='__main__':main()
