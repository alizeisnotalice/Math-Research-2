#!/usr/bin/env python3
"""R59 complete boundary-halfplane surgery: endpoint mass split, no N scan."""
import argparse,bisect,hashlib,json,math,os,platform,sys,time
from fractions import Fraction as F
from pathlib import Path
from periodic_halo import union_area


def ceilq(x):return -(-x.numerator//x.denominator)
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def rootup(x,den):
    k=math.isqrt(x.numerator*den*den//x.denominator)
    assert F(k,den)**2<=x<F(k+1,den)**2
    return F(k+1,den)


def setup(data,side,t):
    px,py,lam=(F(data[k]) for k in ('Px','Py','lam'))
    raw=[tuple(map(F,p)) for p in data['pts']]
    assert all(0<=x<px and 0<=y<py and m>0 for x,y,m in raw)
    assert 0<=t<=1 and lam==1
    X=max(x for x,y,m in raw);labels=[j for j,(x,y,m) in enumerate(raw) if x==0]
    assert labels
    if side=='left':
        pts=[(X-x,y,m,t if j in labels else F(1)) for j,(x,y,m) in enumerate(raw)]
        anchor=X-3*px
    elif side=='right':
        pts=[(x,y,m,1-t if j in labels else F(0)) for j,(x,y,m) in enumerate(raw)]
        anchor=-3*px
    else:raise ValueError(side)
    mass=sum(m for x,y,m in raw);rho=mass/(px*py)
    assert rho*(1+px/20)*(1+py/20)<1
    assert rho*(1+px/15)*(1+py/30)<1
    xmax=max(x if z>0 else x-px for x,y,m,z in pts)
    return px,py,lam,pts,anchor,xmax,{'side':side,'t':str(t),'selected_labels':labels,'selected_mass':str(sum(raw[j][2] for j in labels)),'X':str(X),'cell_mass':str(mass),'density':str(rho),'anchor':str(anchor),'actual_source_xmax':str(xmax)}


def source_atoms(pts,px,py,anchor,C,R):
    atoms={};tags=[]
    for mirror in [False,True]:
        lx,ux=(C-R,2*C-anchor+R) if mirror else (anchor-R,C+R)
        for label,(a,b,m,z) in enumerate(pts):
            imin=ceilq((lx-a)/px);imax=min(0,(ux-a)//px)
            jmin=ceilq((-R-b)/py);jmax=(py+R-b)//py
            for i in range(imin,imax+1):
                mass=m*(z if i==0 else 1)
                if mass==0:continue
                x=a+i*px
                if mirror:x=2*C-x
                for j in range(jmin,jmax+1):
                    y=b+j*py;atoms[x,y]=atoms.get((x,y),F(0))+mass
                    tags.append((mirror,label,i,j))
    return atoms,tags


def prefix(atoms,scale,W):
    raw={(int(x*scale),int(y*scale)):int(m*W) for (x,y),m in atoms.items()}
    xs=sorted({x for x,y in raw});ys=sorted({y for x,y in raw});nx,ny=len(xs),len(ys)
    xm={x:i+1 for i,x in enumerate(xs)};ym={y:j+1 for j,y in enumerate(ys)}
    p=[[0]*(ny+1) for _ in range(nx+1)]
    for (x,y),m in raw.items():p[xm[x]][ym[y]]+=m
    for i in range(1,nx+1):
        acc=0
        for j in range(1,ny+1):acc+=p[i][j];p[i][j]=p[i-1][j]+acc
    return raw,xs,ys,p


def axis(data,side,t,bits):
    px,py,lam,pts,anchor,xmax,meta=setup(data,side,t)
    R=F(10);boot=[]
    for step in range(40):
        nr=ceilq(2*R/py)+1;row=[]
        for label,(a,b,m,z) in enumerate(pts):
            imin=(xmax-R-a)//px+1
            for i in range(imin,1):
                mm=m*(z if i==0 else 1)
                if mm>0:row.append((i,label,a+i*px,b,mm))
        B=nr*sum(p[4] for p in row);new=rootup(B/4,1000)
        boot.append({'R':str(R),'rows':nr,'mass_upper':str(B),'next_strict_radius':str(new),'row_labels':[[i,label,str(m)] for i,label,x,y,m in row]})
        if new>=R:break
        R=new
    else:raise AssertionError('bootstrap did not stop')
    atoms={}
    for i,label,x,y,m in row:
        for j in range(nr):atoms[x,y+j*py]=atoms.get((x,y+j*py),F(0))+m
    D=1<<bits;S=math.lcm(D,R.denominator,*(z.denominator for p in atoms for z in p));W=math.lcm(*(m.denominator for m in atoms.values()))
    raw,xs,ys,p=prefix(atoms,S,W);Ri=int(R*S);B=4*W;best=None;owner=None
    count={'boxes':0,'span_rejected':0,'cannot_beat_lower':0,'root_checked':0}
    for (x,y),m in raw.items():
        rr=math.isqrt(m*D*D//B);v=x+rr*(S//D)
        if best is None or v>best:best=v;owner={'bbox':[str(F(x,S)),str(F(x,S)),str(F(y,S)),str(F(y,S))],'mass':str(F(m,W))}
    for i,a in enumerate(xs):
      for k in range(i,bisect.bisect_left(xs,a+2*Ri)):
        b=xs[k]
        for j,c in enumerate(ys):
          for l in range(j,bisect.bisect_left(ys,c+2*Ri)):
            d=ys[l];count['boxes']+=1;M=p[k+1][l+1]-p[i][l+1]-p[k+1][j]+p[i][j]
            if not M or max(b-a,d-c)**2*B>=4*S*S*M:count['span_rejected']+=1;continue
            assert best>a
            if S*S*M<=B*(best-a)**2:count['cannot_beat_lower']+=1;continue
            rr=math.isqrt(M*D*D//B);v=a+rr*(S//D);count['root_checked']+=1
            if v>best:best=v;owner={'bbox':[str(F(a,S)),str(F(b,S)),str(F(c,S)),str(F(d,S))],'mass':str(F(M,W))}
    lo,hi=F(best,S),F(best,S)+F(1,D)
    assert lo>xmax
    return {'status':'PASS','mode':'axis',**meta,'bootstrap':boot,'final_owner_radius':str(R),'catalog_atoms':len(atoms),'catalog_x':len(xs),'catalog_y':len(ys),'counts':count,'c_lower':str(lo),'c_upper':str(hi),'rational_axis':str(hi),'delta':str(hi-lo),'max_lower_witness':owner,'bits':bits,'scope':'Complete right halfplane source-bbox maximum enclosure, with labelwise modified boundary column; left side uses exact nonmodular xprime=X-x.'}


def area(data,side,t,seam,bits):
    px,py,lam,pts,anchor,xmax,meta=setup(data,side,t)
    assert seam['mode']=='axis' and seam['side']==side and F(seam['t'])==t
    lo,hi=F(seam['c_lower']),F(seam['c_upper']);C=F(seam['rational_axis']);assert C==hi and lo>xmax
    R=F(15);rho=F(meta['density']);tail=rho*(1+px/15)*(1+py/30);assert tail<1
    atoms,tags=source_atoms(pts,px,py,anchor,C,R)
    D=1<<bits;scale=math.lcm(D,anchor.denominator,C.denominator,py.denominator,*(z.denominator for p in atoms for z in p));W=math.lcm(*(m.denominator for m in atoms.values()))
    raw,xs,ys,p=prefix(atoms,scale,W);nx,ny=len(xs),len(ys);B=4*W;Ri=int(R*scale)
    X0,X1,Y0,Y1=int(anchor*scale),int(C*scale),0,int(py*scale)
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
            count['boxes']+=1;M=rowhi[l+1]-rowlo[l+1]-rowhi[j]+rowlo[j]
            if not M:count['zero_mass']+=1;continue
            span=max(b-a,d-c)
            if span*span*B>=4*scale2*M:count['span_rejected']+=1;continue
            rr=math.isqrt(M*D2//B);count['root_checked']+=1
            assert rr*rr*B<=M*D2<(rr+1)*(rr+1)*B
            for radius,target,label in [(min(Ri,max(0,rr-1)*factor),inner,'inner_empty'),(min(Ri,(rr+1)*factor),outer,'outer_empty')]:
                q=max(X0,b-radius),min(X1,a+radius),max(Y0,d-radius),min(Y1,c+radius)
                if q[0]<q[1] and q[2]<q[3]:target.add(q)
                else:count[label]+=1
    ai,ao=F(union_area(inner),scale2),F(union_area(outer),scale2)
    common,ctags=source_atoms(pts,px,py,anchor,lo,R)
    assert set(tags)<=set(ctags)
    cnx=len({x for x,y in common});cny=len({y for x,y in common});allbox=cnx*(cnx+1)//2*cny*(cny+1)//2
    error=2*(hi-lo)*py*allbox
    return {'status':'PASS','mode':'area',**meta,'window':[str(anchor),str(C),'0',str(py)],'tail_radius':'15','tail_bound':str(tail),'rational_axis':str(C),'bits':bits,'scale':str(scale),'weight_scale':str(W),'window_atoms':len(atoms),'window_x':nx,'window_y':ny,'counts':count,'inner_area':str(ai),'outer_area':str(ao),'inner_area_float_display':float(ai),'outer_area_float_display':float(ao),'inner_rectangles':len(inner),'outer_rectangles':len(outer),'true_axis_error':str(error),'true_axis_area_lower':str(ai-error),'true_axis_area_upper':str(ao+error),'common_catalog':{'atoms':len(common),'nx':cnx,'ny':cny,'all_bbox_upper_count':allbox,'high_axis_tags_subset_low_axis_tags':True},'scope':'Full reflected-halfplane halo in anchored original-side collar; complete source bboxes, cap15, exact union, all-bbox true-axis perturbation including moving seam window.'}


def main():
    ap=argparse.ArgumentParser();ap.add_argument('input');ap.add_argument('--side',choices=['left','right'],required=True);ap.add_argument('--t',default='1/2');ap.add_argument('--mode',choices=['axis','area'],required=True);ap.add_argument('--axis');ap.add_argument('--bits',type=int);ap.add_argument('--output',required=True);a=ap.parse_args()
    assert os.path.ismount('/Volumes/data');start=time.time();path=Path(a.input);data=json.loads(path.read_text());t=F(a.t)
    if a.mode=='axis':r=axis(data,a.side,t,a.bits or 64)
    else:
        assert a.axis;sp=Path(a.axis);seam=json.loads(sp.read_text());assert seam['input_sha256']==sha(path)
        r=area(data,a.side,t,seam,a.bits or 40);r['axis_sha256']=sha(sp)
    r.update(input_sha256=sha(path),script_sha256=sha(__file__),union_module_sha256=sha(Path(__file__).with_name('periodic_halo.py')),python=sys.version,platform=platform.platform(),elapsed_seconds=time.time()-start,random_seed=None,floating_tolerance=None)
    out=Path(a.output);out.parent.mkdir(parents=True,exist_ok=True);out.write_text(json.dumps(r,indent=2)+'\n')
    print(json.dumps({k:r[k] for k in ['status','mode','side','t','counts','elapsed_seconds']},indent=2))

if __name__=='__main__':main()
