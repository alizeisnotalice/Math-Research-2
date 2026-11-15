#!/usr/bin/env python3
"""Complete fixed-source one-axis halo bbox enumeration, using rational margins.

API: enumerate_slab(data, velocities, s0, H, R, domain=None, periodic=True).
All input numbers should be exact rational strings. Output boxes are raw (unclipped)
(xlo,xhi,ylo,yhi), each endpoint (a,v,sign,q) meaning a+v*(s-s0)+sign*sqrt(q).
Only source/span/reachability-stable bboxes are prefiltered. Every subsequent
empty bbox has a permanent reason or a positive rational exclusion radius.
No endpoint ordering or slope computation is claimed here.
"""
import argparse
import hashlib
import json
import math
import os
import platform
import sys
import time
from collections import Counter
from fractions import Fraction as F
from pathlib import Path


class WallAtCenter(ValueError):
    """A nonpersistent event at s0 prevents a symmetric stability certificate."""
    def __init__(self, kind, witness):
        self.kind, self.witness = kind, witness
        super().__init__(f'{kind}: {witness}')


def _ceil(q):
    return -((-q.numerator)//q.denominator)


def _sha(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def _expr(a, v=0, sign=0, q=0):
    return [str(a), str(v), int(sign), str(q)]


def _radius_expr(a, v, sign, r2):
    n, d = math.isqrt(r2.numerator), math.isqrt(r2.denominator)
    if n*n == r2.numerator and d*d == r2.denominator:
        return _expr(a + sign*F(n,d), v)
    return _expr(a, v, sign, r2)


def _serial(obj):
    if isinstance(obj,F): return str(obj)
    if isinstance(obj,dict): return {k:_serial(v) for k,v in obj.items()}
    if isinstance(obj,(list,tuple)): return [_serial(v) for v in obj]
    return obj


def enumerate_slab(data, velocities, s0, H, R, *, domain=None, periodic=True):
    """Return a complete local enumeration or raise WallAtCenter / ValueError.

    H initially bounds the entire fixed periodic copy catalogue. It is shrunk to
    half the distance to every nonpersistent source/order/span/reach event. The
    output delta_empty further controls all empty rectangles. No omitted bbox
    can become nonempty within delta_empty. Endpoint sweep must additionally
    certify all active raw endpoint order/clip branches on its final radius.
    """
    start = time.perf_counter()
    s0,H,R,lam = map(F,(s0,H,R,data['lam']))
    if min(H,R,lam)<=0: raise ValueError('H,R,lambda must be positive')
    pts = [tuple(map(F,p)) for p in data['pts']]
    if len(velocities)!=len(pts) or any(F(v) not in (0,1) for v in velocities):
        raise ValueError('one 0/1 velocity per original atom required')
    velocities=list(map(int,map(F,velocities)))
    if any(w<0 for _,_,w in pts) or not sum(p[2] for p in pts):
        raise ValueError('nonnegative weights and positive total mass required')
    px,py=(F(data[k]) for k in ('Px','Py')) if periodic else (None,None)
    if periodic and min(px,py)<=0: raise ValueError('positive periods required')
    domain=tuple(map(F,domain or (0,px,0,py)))
    fl,fr,fb,ft=domain
    if fl>=fr or fb>=ft: raise ValueError('positive domain required')
    mass=sum(p[2] for p in pts)
    tail=mass/(px*py)*(1+px/(2*R))*(1+py/(2*R)) if periodic else mass/(4*R*R)
    if tail>=lam: raise ValueError('position-independent tail bound must be strictly below lambda')
    initial_H=H
    sources=[]
    for label,((x,y,w),v) in enumerate(zip(pts,velocities)):
        if not w: continue
        xc=x+v*s0
        if periodic:
            xmin,xmax=xc-v*H,xc+v*H
            ir=range(_ceil((fl-R-xmax)/px), (fr+R-xmin)//px+1)
            jr=range(_ceil((fb-R-y)/py), (ft+R-y)//py+1)
        else: ir=jr=(0,)
        for i in ir:
            xx=xc+(i*px if periodic else 0)
            for j in jr:
                yy=y+(j*py if periodic else 0)
                sources.append({'label':label,'period_shift':[i,j], 'x':xx,'y':yy,'v':v,'w':w})
    if not sources: raise ValueError('empty source catalogue')
    # Coordinate classes combine only persistent affine coordinates. Distinct
    # velocities meeting at the center are explicitly rejected as source walls.
    xclasses=sorted({(p['x'],p['v']) for p in sources})
    ys=sorted({p['y'] for p in sources})
    source_limiter={'kind':'requested_slab','radius':H}
    event_counts=Counter()
    def event(gap, speed, kind, witness):
        nonlocal H,source_limiter
        if not speed:
            event_counts[kind+'_persistent_equal' if not gap else kind+'_parallel']+=1
            return
        if not gap: raise WallAtCenter(kind, _serial(witness))
        radius=abs(gap/speed)/2
        event_counts[kind+'_nonpersistent']+=1
        if radius<H:
            H=radius;source_limiter={'kind':kind,'radius':H,'gap':gap,'relative_speed':speed,'witness':witness}
    for i,(a,va) in enumerate(xclasses):
        for k in range(i+1,len(xclasses)):
            b,vb=xclasses[k]
            event(b-a,vb-va,'source_order',{'xclasses':[i,k]})
            event(b-a-2*R,vb-va,'span_2R',{'xclasses':[i,k]})
        # These are all inequalities actually used for x reach prefiltering.
        event(a-(fr+R),va,'reach_right',{'xclass':i,'source_x':a})
        event(a-(fl-R),va,'reach_left',{'xclass':i,'source_x':a})
    if H<=0: raise ValueError('no positive source-order slab')
    # Persistent class coordinates are strictly increasing after collision check.
    xs=[z[0] for z in xclasses]; vs=[z[1] for z in xclasses]
    if any(a>=b for a,b in zip(xs,xs[1:])): raise ValueError('nonpersistent collision escaped event test')
    nx,ny=len(xs),len(ys)
    scale=math.lcm(*(z.denominator for z in [*xs,*ys,*domain,R]))
    weight_scale=math.lcm(*(p['w'].denominator for p in sources))
    X=[int(z*scale) for z in xs];Y=[int(z*scale) for z in ys]
    RI=int(R*scale); FL,FR,FB,FT=(int(z*scale) for z in domain)
    ln,ld=lam.numerator,lam.denominator
    xi={p:i for i,p in enumerate(xclasses)};yi={p:i for i,p in enumerate(ys)}
    prefix=[[0]*(ny+1) for _ in range(nx+1)]
    for p in sources:prefix[xi[p['x'],p['v']]+1][yi[p['y']]+1]+=int(p['w']*weight_scale)
    for i in range(1,nx+1):
        acc=0
        for j in range(1,ny+1):
            acc+=prefix[i][j];prefix[i][j]=prefix[i-1][j]+acc
    xpairs=[];ypairs=[];counts=Counter()
    for i,a in enumerate(X):
        for k in range(i,nx):
            b=X[k]
            if b-a>2*RI:counts['x_pair_span_prefilter']+=1;continue
            if b-RI>=FR or a+RI<=FL:counts['x_pair_reach_prefilter']+=1;continue
            xpairs.append((i,k))
    for j,c in enumerate(Y):
        for l in range(j,ny):
            d=Y[l]
            if d-c>2*RI:counts['y_pair_span_prefilter']+=1;continue
            if d-RI>=FT or c+RI<=FB:counts['y_pair_reach_prefilter']+=1;continue
            ypairs.append((j,l))
    hn,hd=H.numerator,H.denominator
    delta=H;minimum={'kind':'source_slab','radius':H,'witness':source_limiter}
    witnesses={};boxes=[];records=[]
    def witness(ids,M,margin=None):
        i,k,j,l=ids
        w={'bbox':list(ids),'mass':F(M,weight_scale),'x_endpoints':[xs[i],xs[k]],
           'x_velocities':[vs[i],vs[k]],'y_endpoints':[ys[j],ys[l]]}
        if margin:
            mode,Q,u,D,B=margin
            w['q']=F(Q,ld*scale*scale*weight_scale) if mode=='mass' else F(Q,scale*scale*B)
            w['L']=(lam if mode=='mass' else 1)*abs(u)*(2*abs(F(D,scale))+abs(u)*H)
            w['margin_kind']=mode
        return w
    def reject(kind,ids,M,radius=None,margin=None):
        nonlocal delta,minimum
        counts[kind]+=1
        if radius is None:
            if kind not in witnesses:witnesses[kind]={'permanent':True,'witness':witness(ids,M)}
            return
        rn,rd=radius
        if rn<=0 or rd<=0:raise ValueError('nonpositive exclusion radius')
        if rn*hd>hn*rd:rn,rd=hn,hd
        old=witnesses.get(kind)
        if old is None or rn*old['radius'].denominator<old['radius'].numerator*rd:
            w=witness(ids,M,margin);rad=F(rn,rd);witnesses[kind]={'radius':rad,'witness':w}
            if rn*delta.denominator<delta.numerator*rd:
                delta=rad;minimum={'kind':kind,'radius':rad,'witness':w}
    # All inner-loop positivity tests and prefix masses use Python integers.
    # Fractions and endpoint square roots are constructed only for survivors
    # or a newly strongest rejection witness, not for every rejected bbox.
    B=4*weight_scale*ln
    for i,k in xpairs:
      a,b=X[i],X[k];va,vb=vs[i],vs[k];D=b-a;u=vb-va
      rowhi,rowlo=prefix[k+1],prefix[i]
      for j,l in ypairs:
        counts['boxes_visited']+=1
        ids=(i,k,j,l);c,d=Y[j],Y[l]
        M=rowhi[l+1]-rowlo[l+1]-rowhi[j]+rowlo[j]
        if not M:reject('zero_mass',ids,M);continue
        # Effective radius in integer coordinate units: sqrt(K/B).
        K=min(M*ld*scale*scale,RI*RI*B)
        if (d-c)**2*B>=4*K:reject('y_span_permanent',ids,M);continue
        if (c<=FB and (FB-c)**2*B>=K) or (d>=FT and (d-FT)**2*B>=K):
            reject('y_clip_permanent',ids,M);continue
        Q=ln*D*D*weight_scale-M*ld*scale*scale
        if Q>=0:
            if Q==0:
                if u:raise WallAtCenter('x_mass_span',_serial(witness(ids,M)))
                reject('x_mass_span_persistent_equal',ids,M);continue
            if not u:reject('x_mass_span_permanent',ids,M);continue
            # q/(2L), with every denominator cleared exactly.
            rad=(Q*hd,2*weight_scale*ln*scale*abs(u)*(2*abs(D)*hd+abs(u)*hn*scale))
            reject('x_mass_span',ids,M,rad,('mass',Q,u,D,B));continue
        if D*D*B>=4*K:
            if u:raise WallAtCenter('x_cap_span',_serial(witness(ids,M)))
            reject('x_cap_span_permanent',ids,M);continue
        clipped=False
        for dist,speed,kind in ((FL-a,-va,'x_clip_left'),(b-FR,vb,'x_clip_right')):
            if dist<0:continue
            Qc=dist*dist*B-K
            if Qc<0:continue
            if Qc==0:
                if speed:raise WallAtCenter(kind,_serial(witness(ids,M)))
                reject(kind+'_persistent_equal',ids,M)
            elif not speed:reject(kind+'_permanent',ids,M)
            else:
                rad=(Qc*hd,2*B*scale*abs(speed)*(2*abs(dist)*hd+abs(speed)*hn*scale))
                reject(kind,ids,M,rad,('clip_square',Qc,speed,dist,B))
            clipped=True;break
        if clipped:continue
        r2=F(K,B*scale*scale)
        boxes.append([_radius_expr(xs[k],vb,-1,r2),_radius_expr(xs[i],va,1,r2),
                      _radius_expr(ys[l],0,-1,r2),_radius_expr(ys[j],0,1,r2)])
        records.append({'bbox':list(ids),'mass':str(F(M,weight_scale)),
                        'radius_squared':str(r2),'x_velocities':[va,vb]})
        counts['active_boxes']+=1
    if counts['boxes_visited']!=len(xpairs)*len(ypairs):raise AssertionError('incomplete bbox loop')
    rejected=sum(v for k,v in counts.items() if k.startswith(('zero_mass','y_span','y_clip','x_mass','x_cap','x_clip')))
    if counts['active_boxes']+rejected!=counts['boxes_visited']:raise AssertionError('unaccounted boxes')
    out={'schema':'position-slab-enumerator-v1','status':'complete','scope':'fixed catalogue full coordinate-bbox cover; active endpoint order still requires separate certification',
         'periodic':periodic,'domain':domain,'s0':s0,'requested_H':initial_H,'H':H,'delta_empty':delta,
         'tail':{'radius':R,'bound':tail,'lambda':lam,'mass':mass,'strict':True,'method':'periodic_count' if periodic else 'finite_total_mass'},
         'input_canonical_sha256':_sha(data),'direction_sha256':_sha(velocities),'velocities':velocities,
         'source_catalog':sources,'source_catalog_sha256':_sha(_serial(sources)),
         'x_classes':xclasses,'y_classes':ys,'source_order_sha256':_sha(_serial([xclasses,ys])),
         'source_event_counts':dict(event_counts),'source_limiter':source_limiter,
         'counts':dict(counts),'candidate_x_pairs':len(xpairs),'candidate_y_pairs':len(ypairs),
         'coordinate_scale':str(scale),'weight_scale':str(weight_scale),'arithmetic':'integer bbox loop; rational event/margin certificates',
         'all_bboxes_before_prefilter':nx*(nx+1)//2*ny*(ny+1)//2,
         'rejection_witnesses':witnesses,'min_empty_witness':minimum,'boxes':boxes,'box_records':records,
         'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
         'python':sys.version,'platform':platform.platform(),'elapsed_seconds':time.perf_counter()-start}
    return _serial(out)


def _oracle_selftest():
    # Independent receiver-point criterion: order atom distances and test the
    # cumulative mass against 4*lambda*r^2, without using bbox geometry.
    import random
    rng=random.Random(39015);checks=0;cases=[]
    for case in range(12):
        pts=[[str(F(rng.randrange(-20,21),7)),str(F(rng.randrange(-20,21),7)),str(rng.randrange(1,4))] for _ in range(5)]
        velocities=[rng.randrange(2) for _ in pts];s0=F(case+1,10007)
        out=enumerate_slab({'lam':'1','pts':pts},velocities,s0,F(1,20000),3,domain=['-4','4','-4','4'],periodic=False)
        dt=F(out['delta_empty']);xs=[(F(x),int(v)) for x,v in out['x_classes']];ys=list(map(F,out['y_classes']))
        for t in [-dt,F(0),dt]:
            pp=[(F(x)+v*(s0+t),F(y),F(w)) for (x,y,w),v in zip(pts,velocities)]
            for _ in range(35):
                x,y=F(rng.randrange(-399,400),100),F(rng.randrange(-399,400),100)
                distances=sorted((max(abs(xx-x),abs(yy-y)),w) for xx,yy,w in pp)
                direct=False;cum=F(0)
                for r,w in distances:
                    cum+=w
                    if r<3 and cum>4*r*r:direct=True
                covered=False
                for rec in out['box_records']:
                    i,k,j,l=rec['bbox'];r2=F(rec['radius_squared'])
                    a=xs[i][0]+xs[i][1]*t;b=xs[k][0]+xs[k][1]*t
                    dist=max(abs(x-a),abs(x-b),abs(y-ys[j]),abs(y-ys[l]))
                    if dist*dist<r2:covered=True;break
                assert direct==covered,(case,t,x,y,direct,covered)
                checks+=1
        cases.append({'pts':pts,'velocities':velocities,'s0':str(s0),'delta_empty':out['delta_empty']})
    return {'seed':39015,'checks':checks,'cases':cases,'criterion':'exact rational cumulative mass ordered by receiver distance'}


def small_selftest():
    """Actual exact checks: collision, death, clipping, cap, catalogue, margins."""
    data={'lam':'1','pts':[['0','0','9'],['0','0','16']]}
    domain=['-10','11','-4','4'];results=[]
    for s in ['-1/2','1/2','2','7/2','17/4','9/2','49/10','51/10','13/2']:
        out=enumerate_slab(data,[0,1],s,'1/10','3',domain=domain,periodic=False)
        # For two atoms exactly 3 bboxes exist away from collision. At s>5 the
        # joint rectangle is rational-margin excluded, with positive radius.
        assert out['counts']['boxes_visited']==(3 if abs(F(s))<6 else 2)
        expected=2 if abs(F(s))>5 else 3
        assert len(out['boxes'])==expected
        assert F(out['delta_empty'])>0
        # Active boxes may hit walls within delta_empty: their endpoint order
        # is deliberately the separate sweep's responsibility.
        if abs(F(s))>5:
            for t in [-F(out['delta_empty']),F(0),F(out['delta_empty'])]:
                assert abs(F(s)+t)>5
        results.append({'s0':s,'H':out['H'],'delta_empty':out['delta_empty'],'active_boxes':len(out['boxes'])})
    for s,kind in [('0','source_order'),('5','x_mass_span'),('6','span_2R')]:
        try:enumerate_slab(data,[0,1],s,'1/10','3',domain=domain,periodic=False)
        except WallAtCenter as e: assert e.kind==kind
        else:raise AssertionError('wall did not fail closed')
    # Clip rejection can be strict at s0 with a nonzero velocity; prove it
    # remains empty without opening the radical of this bbox.
    clip={'lam':'1','pts':[['0','0','9'],['3','0','1']]}
    out=enumerate_slab(clip,[0,1],0,'1/10',3,domain=['-1','1','-2','2'],periodic=False)
    assert len(out['boxes'])==1 and out['counts']['x_clip_right']>=1
    assert F(out['delta_empty'])>0
    clip['pts'][1][0]='3/2'
    try:enumerate_slab(clip,[0,1],0,'1/10',3,domain=['-1','1','-2','2'],periodic=False)
    except WallAtCenter as e:assert e.kind=='x_clip_right'
    else:raise AssertionError('clip wall did not fail closed')
    # A periodic catalogue contains every copy needed at both initial endpoints.
    periodic={'Px':'20','Py':'20','lam':'1','pts':[['0','0','1'],['19','19','1']]}
    out=enumerate_slab(periodic,[1,0],'1/20','1/10',3)
    catalog={(p['label'],*p['period_shift']) for p in out['source_catalog']}
    for t in [-F(1,10),F(1,10)]:
        for label,((xx,yy,_),v) in enumerate(zip(periodic['pts'],[1,0])):
            for i in range(-2,3):
                for j in range(-2,3):
                    x,y=F(xx)+v*(F(1,20)+t)+i*20,F(yy)+j*20
                    if -3<=x<=23 and -3<=y<=23: assert (label,i,j) in catalog
    return {'status':'passed','two_atom_samples':results,'walls_rejected':['source_collision_0','joint_death_5','span_cap_6'],
            'periodic_catalog_endpoint_checks':True,'independent_oracle':_oracle_selftest(),
            'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'python':sys.version,'scope':'finite exact enumeration tests; endpoint slope handled by separate sweep'}


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('input',nargs='?');ap.add_argument('--velocities',help='JSON file containing original-label 0/1 array')
    ap.add_argument('--s0',default='0');ap.add_argument('--H',default='1/100');ap.add_argument('--R',default='10')
    ap.add_argument('--domain',nargs=4);ap.add_argument('--finite',action='store_true');ap.add_argument('--self-test',action='store_true')
    ap.add_argument('--output');args=ap.parse_args()
    if args.self_test:out=small_selftest()
    else:
        if not args.input or not args.velocities:ap.error('input and --velocities required')
        out=enumerate_slab(json.loads(Path(args.input).read_text()),json.loads(Path(args.velocities).read_text()),args.s0,args.H,args.R,domain=args.domain,periodic=not args.finite)
        out['input_file_sha256']=hashlib.sha256(Path(args.input).read_bytes()).hexdigest()
        out['direction_file_sha256']=hashlib.sha256(Path(args.velocities).read_bytes()).hexdigest()
    if args.output:
        if not os.path.ismount('/Volumes/data'):raise RuntimeError('/Volumes/data is not mounted')
        target=Path(args.output);target.parent.mkdir(parents=True,exist_ok=True);target.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out if args.self_test else {k:v for k,v in out.items() if k not in ('boxes','box_records','source_catalog')},indent=2))


if __name__=='__main__':main()
