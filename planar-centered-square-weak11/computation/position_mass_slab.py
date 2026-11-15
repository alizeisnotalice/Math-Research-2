#!/usr/bin/env python3
"""Complete R52 fixed-source position/mass bbox rescan with rational empty margins.

The center data already contains center masses. Its x coordinates are the base
coordinates before applying source_slab.s0. Mass velocities are d(weight)/dt;
all returned endpoint expressions use local dz=z-z0 and dt=t-t0:
    [a,v,sign,q0,q1] = a+v*dz+sign*sqrt(q0+q1*dt).
A square q0 is NEVER folded into a if q1 is nonzero. All eligible bboxes are
rescanned with two exact integer prefix sums. The accepted one-axis source slab
provides fixed copy I, source order and all actual geometric prefilter events;
its old mass-fixed empty radius is NOT reused.
"""
import argparse
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import platform
import sys
import time
from collections import Counter
from fractions import Fraction as F

FROZEN_ENUMERATOR_SHA256='acede2c710655d0c179e7bea3ee9483271caaf37968227e7e9fce6442907056b'


def _serial(x):
    if isinstance(x,F):return str(x)
    if isinstance(x,dict):return {k:_serial(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):return [_serial(v) for v in x]
    return x


def _digest(x):
    return hashlib.sha256(json.dumps(_serial(x),sort_keys=True,separators=(',',':')).encode()).hexdigest()


def _expr(a,v,sign,q0=0,q1=0):return [str(a),str(v),int(sign),str(q0),str(q1)]


def build_position_mass_slab(center_data,mass_velocity,x_velocity,source_slab,*,T='1/1000',t0='1/2'):
    start=time.perf_counter();T,t0=F(T),F(t0)
    if T<=0:raise ValueError('positive T required')
    if source_slab.get('schema')!='position-slab-enumerator-v1' or source_slab.get('status')!='complete':
        raise ValueError('requires accepted original one-axis source slab')
    if source_slab.get('code_sha256')!=FROZEN_ENUMERATOR_SHA256:
        raise ValueError('source slab is not the pinned R39 enumerator output')
    if _digest(center_data)!=source_slab['input_canonical_sha256']:
        raise ValueError('center data does not match source slab input')
    pts=[tuple(map(F,p)) for p in center_data['pts']]
    mu=list(map(F,mass_velocity));vx=list(map(F,x_velocity))
    if len(mu)!=len(pts) or len(vx)!=len(pts) or any(v not in (0,1) for v in vx):
        raise ValueError('one mass derivative and 0/1 x velocity per original atom required')
    vx=list(map(int,vx))
    if _digest(vx)!=source_slab['direction_sha256'] or vx!=source_slab['velocities']:
        raise ValueError('x velocity differs from accepted source slab')
    if any(w<=0 for _,_,w in pts):raise ValueError('every center atom mass must be strictly positive')
    if any(w-abs(m)*T<0 for (_,_,w),m in zip(pts,mu)):
        raise ValueError('an atom mass becomes negative on full requested T interval')
    z0,H,R,lam=F(source_slab['s0']),F(source_slab['H']),F(source_slab['tail']['radius']),F(center_data['lam'])
    if min(H,R,lam)<=0:raise ValueError('positive source H/R/lambda required')
    periodic=bool(source_slab['periodic']);domain=tuple(map(F,source_slab['domain']));fl,fr,fb,ft=domain
    px,py=(F(center_data[k]) for k in ('Px','Py')) if periodic else (None,None)
    Mtotal=sum(p[2] for p in pts);MUtotal=sum(mu);Mmax=Mtotal+abs(MUtotal)*T
    tail=Mmax/(px*py)*(1+px/(2*R))*(1+py/(2*R)) if periodic else Mmax/(4*R*R)
    if tail>=lam:raise ValueError('uniform full-T total-mass tail bound is not strict')
    # This ensures every included root stays positive even if an atom vanishes
    # at a full-T endpoint; no zero derivative/square-root division is assumed.
    rho_mass=F(1);mass_limiter=None
    for label,((_,_,w),m) in enumerate(zip(pts,mu)):
        if m:
            rr=min(F(1),w/(2*abs(m)*T))
            if rr<rho_mass:rho_mass=rr;mass_limiter={'label':label,'w0':w,'mu':m,'rho':rr}
    sources=source_slab['source_catalog']
    if _digest(sources)!=source_slab['source_catalog_sha256']:raise ValueError('source catalogue hash mismatch')
    for source in sources:
        label=source['label'];i,j=source['period_shift'];x,y,w=pts[label]
        expected_x=x+vx[label]*z0+(i*px if periodic else 0)
        expected_y=y+(j*py if periodic else 0)
        if (F(source['x']),F(source['y']),F(source['w']),source['v'])!=(expected_x,expected_y,w,vx[label]):
            raise ValueError('source identity, shift, center geometry or mass binding mismatch')
    classes=[(F(x),int(v)) for x,v in source_slab['x_classes']];ys=list(map(F,source_slab['y_classes']))
    if _digest([source_slab['x_classes'],source_slab['y_classes']])!=source_slab['source_order_sha256']:
        raise ValueError('source order signature mismatch')
    if any(a[0]>=b[0] for a,b in zip(classes,classes[1:])) or any(a>=b for a,b in zip(ys,ys[1:])):
        raise ValueError('source coordinate class order not strict')
    nx,ny=len(classes),len(ys);xs=[x for x,v in classes];vs=[v for x,v in classes]
    S=math.lcm(*(q.denominator for q in [*xs,*ys,*domain,R]))
    W=math.lcm(*(q.denominator for q in [*(p[2] for p in pts),*mu]))
    X=[int(x*S) for x in xs];Y=[int(y*S) for y in ys]
    RI=int(R*S);FL,FR,FB,FT=(int(q*S) for q in domain)
    ix={q:i for i,q in enumerate(classes)};iy={q:i for i,q in enumerate(ys)}
    p0=[[0]*(ny+1) for _ in range(nx+1)];p1=[[0]*(ny+1) for _ in range(nx+1)]
    for source in sources:
        i=ix[F(source['x']),source['v']]+1;j=iy[F(source['y'])]+1;label=source['label']
        p0[i][j]+=int(pts[label][2]*W);p1[i][j]+=int(mu[label]*W)
    for i in range(1,nx+1):
        a=b=0
        for j in range(1,ny+1):
            a+=p0[i][j];b+=p1[i][j];p0[i][j]=p0[i-1][j]+a;p1[i][j]=p1[i-1][j]+b
    counts=Counter();xpairs=[];ypairs=[]
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
    if len(xpairs)*len(ypairs)!=source_slab['counts']['boxes_visited']:
        raise ValueError('eligible bbox count differs from accepted source slab')
    ln,ld=lam.numerator,lam.denominator;B=4*W*ln
    hn,hd=H.numerator,H.denominator;tn,td=T.numerator,T.denominator
    RN=hd*hd*td;AZ=B*S*hn*td;AM=S*S*ld*tn*hd*hd
    empty_num=empty_den=1;empty_limiter={'kind':'requested_rho_one','rho':'1'}
    cap_rho=F(1);cap_limiter=None;witnesses={};boxes=[];records=[];active_keys=[];wall_keys=[]
    def margin(kind,D,u,j,clip,M,MU,K):
        if D<0:return None
        Q=D*D*B-j*j*K
        if Q<0:return None
        if not Q:
            cap_constant=(M*td-abs(MU)*tn)*ld*S*S>=RI*RI*B*td
            if not u and (not MU or cap_constant):return (1,1,kind,D,u,j,Q,True,clip)
            return (0,1,kind,D,u,j,Q,False,clip)
        L=abs(u)*(2*abs(D)*hd+abs(u)*hn*S)*AZ+j*j*abs(MU)*AM
        if not L:return (1,1,kind,D,u,j,Q,True,clip)
        rn,rd=Q*RN,2*L
        if clip and u:
            cn,cd=abs(D)*hd,2*abs(u)*hn*S
            if cn*rd<rn*cd:rn,rd=cn,cd
        if rn>=rd:rn=rd=1
        return (rn,rd,kind,D,u,j,Q,False,clip)
    def witness(ids,M,MU,choice):
        rn,rd,kind,D,u,j,Q,permanent,clip=choice
        d0=F(D,j*S);velocity=F(u,j);r02=min(R*R,F(M*ld,B))
        q=d0*d0-r02;L=abs(velocity)*(2*abs(d0)+abs(velocity)*H)*H+F(abs(MU)*ld,B)*T
        return {'bbox':list(ids),'m0':F(M,W),'mu':F(MU,W),'reason':kind,'d0':d0,'d_velocity':velocity,
                'r0_squared':r02,'q':q,'L':L,'rho_bound':F(rn,rd),'permanent':permanent,
                'clip_distance_positivity_required':clip,'clip_positive_radius':d0/(2*abs(velocity)*H) if clip and velocity else None,
                'guarantee':'q(z,t)>=q0/2>0; clip d remains >=d0/2 when nonconstant' if q>0 else 'persistent equality of the failed strict condition'}
    def reject(ids,M,MU,choice):
        nonlocal empty_num,empty_den,empty_limiter
        rn,rd,kind,*_=choice;counts['reject_'+kind]+=1
        old=witnesses.get(kind)
        if old is None or rn*old[1]<old[0]*rd:
            witnesses[kind]=(rn,rd,witness(ids,M,MU,choice))
        if rn*empty_den<empty_num*rd:
            empty_num,empty_den=rn,rd;empty_limiter=witness(ids,M,MU,choice)
    for i,k in xpairs:
      a,b=X[i],X[k];va,vb=vs[i],vs[k];u=vb-va
      row0,lo0=p0[k+1],p0[i];row1,lo1=p1[k+1],p1[i]
      for j,l in ypairs:
        counts['boxes_visited']+=1;ids=(i,k,j,l);c,d=Y[j],Y[l]
        M=row0[l+1]-lo0[l+1]-row0[j]+lo0[j]
        MU=row1[l+1]-lo1[l+1]-row1[j]+lo1[j]
        if M<0:raise AssertionError('negative center bbox mass')
        if not M:
            if MU:raise ValueError('zero-mass bbox has nonzero mass derivative')
            counts['reject_zero_mass']+=1;continue
        K=min(M*ld*S*S,RI*RI*B)
        best=margin('y_span',d-c,0,2,False,M,MU,K)
        walls=[]
        if best is not None and not best[0]:walls.append(best);best=None
        if best is not None and best[:2]==(1,1):reject(ids,M,MU,best);continue
        for args in (('y_clip_below',FB-c,0,1,True),('y_clip_above',d-FT,0,1,True),
                     ('x_span',b-a,u,2,False),('x_clip_left',FL-a,-va,1,True),('x_clip_right',b-FR,vb,1,True)):
            option=margin(*args,M,MU,K)
            if option is None:continue
            if not option[0]:walls.append(option);continue
            if best is None or option[0]*best[1]>best[0]*option[1]:best=option
            if best[:2]==(1,1):break
        if best is not None:reject(ids,M,MU,best);continue
        C=M*ld*S*S-RI*RI*B
        if not C and MU:raise ValueError('candidate lies on nonpersistent mass-cap branch wall')
        branch='cap' if C>=0 else 'root'
        q0,q1=F(M*ld,B),F(MU*ld,B)
        bound=F(1) if not MU else min(F(1),F(abs(C)*td,2*abs(MU)*ld*S*S*tn))
        if bound<=0:raise ValueError('candidate cap branch has no positive stability radius')
        if bound<cap_rho:
            cap_rho=bound;cap_limiter={'bbox':list(ids),'q0':q0,'q1':q1,'branch':branch,'rho_bound':bound}
        if branch=='cap':
            raw=[_expr(xs[k]-R,vb,0),_expr(xs[i]+R,va,0),_expr(ys[l]-R,0,0),_expr(ys[j]+R,0,0)]
        else:
            raw=[_expr(xs[k],vb,-1,q0,q1),_expr(xs[i],va,1,q0,q1),_expr(ys[l],0,-1,q0,q1),_expr(ys[j],0,1,q0,q1)]
        center_status='empty_equality_wall' if walls else 'active'
        rec={'bbox':list(ids),'m0':F(M,W),'mu':F(MU,W),'radius_branch':branch,
             'q0':q0,'q1':q1,'cap_branch_rho_bound':bound,'x_velocities':[va,vb],
             'center_status':center_status,'unresolved_zero_margin_reasons':[o[2] for o in walls]}
        boxes.append(raw);records.append(rec);counts['candidate_boxes']+=1
        if walls:wall_keys.append(ids);counts['center_empty_wall_candidates']+=1
        else:active_keys.append(ids);counts['center_active_candidates']+=1
    if counts['boxes_visited']!=len(xpairs)*len(ypairs):raise AssertionError('bbox scan incomplete')
    if counts['candidate_boxes']+sum(v for k,v in counts.items() if k.startswith('reject_'))!=counts['boxes_visited']:
        raise AssertionError('bbox rejection/candidate ledger incomplete')
    old_keys={tuple(r['bbox']) for r in source_slab['box_records']};new_keys=set(active_keys)
    if old_keys!=new_keys:raise ValueError('center-active bbox keys disagree with accepted one-axis slab')
    rho_empty=F(empty_num,empty_den)
    constraints=[('empty_box_margins',rho_empty),('candidate_cap_branches',cap_rho),('strict_atom_mass',rho_mass),('requested_rho_one',F(1))]
    limiter,rho=min(constraints,key=lambda q:q[1])
    if rho<=0 or rho>1:raise AssertionError('invalid uniform rho')
    for record in records:
        q0,q1=record['q0'],record['q1'];low,high=q0-abs(q1)*rho*T,q0+abs(q1)*rho*T
        if low<=0:raise AssertionError('candidate root mass not uniformly strictly positive')
        if record['radius_branch']=='root' and high>=R*R:raise AssertionError('root branch cap not strictly separated')
        if record['radius_branch']=='cap' and low<R*R:raise AssertionError('cap branch not stable')
        record['q_full_domain_lower']=low;record['q_full_domain_upper']=high
        record['cap_gap_lower']=(R*R-high) if record['radius_branch']=='root' else (low-R*R)
    out={'schema':'position-mass-slab-v1','status':'complete_candidate_and_empty_certificate',
         'external_premise':'accepted one-axis source catalogue, source order and all geometric prefilter events are valid on H; this module freshly rescans every eligible bbox for mass variation',
         'z0':z0,'t0':t0,'H':H,'T':T,'rho':rho,'rho_empty':rho_empty,'rho_cap':cap_rho,'rho_atom_positive':rho_mass,
         'dz_radius':rho*H,'dt_radius':rho*T,'rho_limiter':limiter,'rho_constraints':constraints,
         'domain':domain,'periodic':periodic,'lambda':lam,'R':R,'candidate_boxes':boxes,'boxes':boxes,'box_records':records,
         'counts':dict(counts),'candidate_x_pairs':len(xpairs),'candidate_y_pairs':len(ypairs),
         'center_active_keys_match_source_slab':True,'center_active_count':len(active_keys),
         'additional_center_empty_wall_keys':wall_keys,'center_active_keys_sha256':_digest(sorted(active_keys)),
         'empty_box_weakest_witness':empty_limiter,'rejection_witnesses':{k:v[2] for k,v in witnesses.items()},
         'cap_branch_weakest_witness':cap_limiter,'strict_mass_weakest_witness':mass_limiter,
         'uniform_mass':{'center_total_mass':Mtotal,'total_derivative':MUtotal,'full_T_mass_upper':Mmax,
                         'center_atoms_strictly_positive':True,'full_T_atom_masses_nonnegative':True,
                         'rho_domain_atom_masses_at_least_half_center':True},
         'uniform_tail':{'radius':R,'bound':tail,'lambda':lam,'strict':True,'full_T':True,
                         'method':'position-independent periodic count with full-T total-mass upper' if periodic else 'finite full-T total mass'},
         'coordinate_scale':str(S),'weight_derivative_common_scale':str(W),'arithmetic':'two integer prefix sums; integer squared-margin comparisons',
         'source_catalog_sha256':source_slab['source_catalog_sha256'],'source_order_sha256':source_slab['source_order_sha256'],
         'source_enumerator_sha256':source_slab['code_sha256'],'source_slab_canonical_sha256':_digest(source_slab),
         'center_data_canonical_sha256':_digest(center_data),'mass_velocity_sha256':_digest(mu),'x_velocity_sha256':_digest(vx),
         'source_mass_derivative_catalog_sha256':_digest([[s['label'],s['period_shift'],mu[s['label']]] for s in sources]),
         'mass_prefix_sha256':_digest(p0),'derivative_prefix_sha256':_digest(p1),
         'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
         'python':sys.version,'platform':platform.platform(),'elapsed_seconds':time.perf_counter()-start}
    return _serial(out)


def _load_enumerator():
    path=Path(__file__).with_name('position_slab_enumerator.py')
    if hashlib.sha256(path.read_bytes()).hexdigest()!=FROZEN_ENUMERATOR_SHA256:raise ValueError('frozen enumerator hash changed')
    spec=importlib.util.spec_from_file_location('position_mass_frozen_enum',path);module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module


def selftest():
    en=_load_enumerator();data={'lam':'1','pts':[['0','0','9'],['0','0','16']]}
    slab=en.enumerate_slab(data,[0,1],'51/10','1/10',3,domain=['-10','11','-4','4'],periodic=False)
    out=build_position_mass_slab(data,[1,-1],[0,1],slab,T='1/1000')
    assert out['center_active_count']==2 and out['counts']['candidate_boxes']==2
    assert F(out['rho'])>0
    # A rational-square q0 is deliberately retained as a moving radical.
    assert any(F(e[3])==F(9,4) and F(e[4])==F(1,4) and e[2] for box in out['boxes'] for e in box)
    rho,H,T=F(out['rho']),F(out['H']),F(out['T'])
    for dz in [-rho*H,F(0),rho*H]:
        for dt in [-rho*T,F(0),rho*T]:
            span=F(51,10)+dz
            assert span*span>25 # joint mass derivative is zero
            assert all(F(p[2])+m*dt>0 for p,m in zip(data['pts'],[1,-1]))
    # A y-span equality previously permanent at fixed mass becomes a candidate
    # wall when its total mass derivative is nonzero. It must not disappear.
    ydata={'lam':'1','pts':[['0','0','9'],['0','5','16']]}
    yslab=en.enumerate_slab(ydata,[0,0],0,'1/10',3,domain=['-4','4','-4','9'],periodic=False)
    wall=build_position_mass_slab(ydata,[1,0],[0,0],yslab,T='1/1000')
    assert wall['counts']['center_empty_wall_candidates']>=1
    assert any(r['center_status']=='empty_equality_wall' and F(r['mu'])==1 for r in wall['box_records'])
    # Empty geometric boxes have exactly zero derivative as well as zero mass.
    diagonal={'lam':'1','pts':[['0','0','1'],['1','1','1']]}
    ds=en.enumerate_slab(diagonal,[0,0],0,'1/10',3,domain=['-4','4','-4','4'],periodic=False)
    diag=build_position_mass_slab(diagonal,[1,-1],[0,0],ds)
    assert diag['counts']['reject_zero_mass']>0
    try:build_position_mass_slab(diagonal,[2000,-2000],[0,0],ds)
    except ValueError:pass
    else:raise AssertionError('negative full-T atom masses accepted')
    return {'status':'passed','scope':'finite exact position/mass bbox examples only; endpoint joint order is a separate certificate',
            'conserving_two_atom':out,'new_mass_activation_wall':wall,'diagonal_zero_mass':diag,
            'negative_full_T_mass_rejected':True,'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('center_data',nargs='?')
    ap.add_argument('--mass-velocity');ap.add_argument('--x-velocity');ap.add_argument('--source-slab');ap.add_argument('--T',default='1/1000');ap.add_argument('--t0',default='1/2')
    ap.add_argument('--self-test',action='store_true');ap.add_argument('--output',required=True);args=ap.parse_args()
    if args.self_test:out=selftest()
    else:
        paths=[args.center_data,args.mass_velocity,args.x_velocity,args.source_slab]
        if any(x is None for x in paths):ap.error('center_data, --mass-velocity, --x-velocity and --source-slab required')
        data,mv,xv,slab=(json.loads(Path(p).read_text()) for p in paths)
        out=build_position_mass_slab(data,mv,xv,slab,T=args.T,t0=args.t0)
        out['input_file_sha256']={name:hashlib.sha256(Path(p).read_bytes()).hexdigest() for name,p in zip(('center_data','mass_velocity','x_velocity','source_slab'),paths)}
        out['command']=[sys.executable,*sys.argv]
    if not os.path.ismount('/Volumes/data'):raise RuntimeError('/Volumes/data is not mounted')
    target=Path(args.output);target.parent.mkdir(parents=True,exist_ok=True);target.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k in ('schema','status','rho','dz_radius','dt_radius','counts','code_sha256','elapsed_seconds')},indent=2))


if __name__=='__main__':main()
