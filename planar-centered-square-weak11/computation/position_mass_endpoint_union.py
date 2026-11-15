#!/usr/bin/env python3
"""R52 exact endpoint-order and independent direct-merge union derivative.

Input: complete position-mass-slab-v1 candidate/empty certificate. This module
never enumerates bboxes. It supplies the missing joint endpoint order, then
computes beta(t)=dA/dz as an exact rational/root sum, beta_t and d_t(d_z J).
Every radical remains identified by (q0,q1); square q0 is folded ONLY if q1=0.
No production endpoint sweep or union tree is imported.
"""
import argparse
from collections import defaultdict
from decimal import Decimal, localcontext
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import gzip
import json
from math import isqrt
import os
from pathlib import Path
import platform
import sys
import time

sys.dont_write_bytecode=True

def serial(x):
    if isinstance(x,F):return str(x)
    if isinstance(x,dict):return {str(k):serial(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):return [serial(v) for v in x]
    return x

def digest(x):return hashlib.sha256(json.dumps(serial(x),sort_keys=True,separators=(',',':')).encode()).hexdigest()

@lru_cache(None)
def root(q,bits):
    if q<0:raise ValueError('negative radical')
    D=1<<bits
    n=isqrt(q.numerator*q.denominator*D*D)//q.denominator
    assert n*n*q.denominator<=q.numerator*D*D<(n+1)*(n+1)*q.denominator
    lo=F(n,D)
    return lo,lo if lo*lo==q else F(n+1,D)

def expr(e):
    a,v,s,q,b=e;a,v,q,b=map(F,(a,v,q,b));s=int(s)
    if s not in (-1,0,1) or q<0:raise ValueError('invalid expression')
    if not s:return a,v,0,F(0),F(0)
    if not b:
        lo,hi=root(q,0)
        # root(...,0) does not recognize nonintegral rational squares.
        n,d=isqrt(q.numerator),isqrt(q.denominator)
        if n*n==q.numerator and d*d==q.denominator:return a+s*F(n,d),v,0,F(0),F(0)
    return a,v,s,q,b

def center_key(e):
    a,v,s,q,b=e
    n,d=isqrt(q.numerator),isqrt(q.denominator)
    if not s or (n*n==q.numerator and d*d==q.denominator):return a+s*F(n,d),0,F(0)
    return a,s,q

def bound(e,bits,dt=F(0)):
    a,v,s,q,b=e
    if not s:return a,a
    lo,hi=root(q+b*dt,bits)
    return (a+lo,a+hi) if s>0 else (a-hi,a-lo)

def ordered(endpoints,bits):
    bycenter=defaultdict(list)
    for e in endpoints:bycenter[center_key(e)].append(e)
    ties=[es for es in bycenter.values() if len(es)>1]
    if ties:
        raise ValueError('nonpersistent center endpoint equality; distinct full 5-tuples: '+str(serial(ties[0])))
    for precision in range(bits,641,80):
        rows=sorted((bound(e,precision)[0],bound(e,precision)[1],e) for e in endpoints)
        if all(a[1]<b[0] for a,b in zip(rows,rows[1:])):return rows,precision
    raise ValueError('distinct endpoint order unresolved at 640 dyadic bits')

def root_sum_eval(constant,terms,bits,dtlo=F(0),dthi=F(0),derivative=False):
    lo=hi=constant if not derivative else F(0)
    for (q,b),c in terms.items():
        qlo,qhi=sorted((q+b*dtlo,q+b*dthi))
        if qlo<=0:raise ValueError('derivative requires uniformly positive radical')
        rlo=root(qlo,bits)[0];rhi=root(qhi,bits)[1]
        if derivative:
            c=c*b/2
            x,y=F(1,rhi),F(1,rlo)
        else:x,y=rlo,rhi
        a,bnd=(c*x,c*y) if c>=0 else (c*y,c*x)
        lo+=a;hi+=bnd
    return lo,hi

def interval_mul(x,y):return min(a*b for a in x for b in y),max(a*b for a in x for b in y)
def interval_sub(x,y):return x[0]-y[1],x[1]-y[0]
def interval_scale(x,c):return (c*x[0],c*x[1]) if c>=0 else (c*x[1],c*x[0])

def certify(slab,bits=160):
    start=time.process_time()
    if slab.get('schema')!='position-mass-slab-v1' or slab.get('status')!='complete_candidate_and_empty_certificate':
        raise ValueError('requires complete joint candidate/empty certificate')
    Z,T=F(slab['dz_radius']),F(slab['dt_radius'])
    if min(Z,T)<=0:raise ValueError('positive base parameter radii required')
    boxes=[tuple(expr(e) for e in box) for box in slab['boxes']]
    if len(boxes)!=len(slab['box_records']) or len(boxes)!=slab['counts']['candidate_boxes']:
        raise ValueError('candidate count mismatch')
    xl,xr,yl,yr=map(F,slab['domain'])
    fixed=lambda x:expr([x,0,0,0,0])
    X0,X1,Y0,Y1=map(fixed,(xl,xr,yl,yr))
    X={X0,X1};Y={Y0,Y1}
    for l,r,a,b in boxes:
        if a[1] or b[1]:raise ValueError('y source motion in z is outside scope')
        X.update((l,r));Y.update((a,b))
    xrows,xbits=ordered(X,bits);yrows,ybits=ordered(Y,bits);bits=max(bits,xbits,ybits)
    rho=F(1);limiter=None;checked=0;rootlower={};lt={}
    for axis,rows in [('x',xrows),('y',yrows)]:
        for _,_,e in rows:
            if e[2]:
                qlo=e[3]-abs(e[4])*T
                if qlo<=0:raise ValueError('base domain radical not strictly positive')
                rlo=root(qlo,bits)[0]
                if rlo<=0:raise ValueError('root lower bound unresolved')
                rootlower[e]=rlo;lt[e]=abs(e[4])/(2*rlo)
            else:lt[e]=F(0)
        for a,b in zip(rows,rows[1:]):
            gap=b[0]-a[1]
            if gap<=0:raise ValueError('center endpoint interval order not strict')
            rate=abs(b[2][1]-a[2][1])*Z+(lt[a[2]]+lt[b[2]])*T
            checked+=1
            if rate:
                candidate=gap/(2*rate)
                if candidate<rho:
                    rho=candidate;limiter={'axis':axis,'left':a[2],'right':b[2],
                        'center_gap_lower':gap,'scaled_rate_bound':rate,'rho_bound':candidate}
    if not 0<rho<=1:raise AssertionError('joint endpoint radius is not positive')
    zr,tr=rho*Z,rho*T
    xranks={row[2]:i for i,row in enumerate(xrows)};yranks={row[2]:i for i,row in enumerate(yrows)}
    events=[[] for _ in xrows];ybounds=[]
    for i,(l,r,a,b) in enumerate(boxes):
        left,right=max(xranks[l],xranks[X0]),min(xranks[r],xranks[X1])
        bottom,top=max(yranks[a],yranks[Y0]),min(yranks[b],yranks[Y1])
        if left>=right or bottom>=top:
            raise ValueError('candidate empty at center; rectangular stable active-cell proof unavailable for bbox '+str(slab['box_records'][i]['bbox']))
        ybounds.append((bottom,top));events[left].append((1,i));events[right].append((-1,i))
    def merge(intervals):
        s=sorted(intervals)
        if not s:return []
        a,b=s[0];out=[]
        for c,d in s[1:]:
            if c>b:out.append((a,b));a,b=c,d
            else:b=max(b,d)
        out.append((a,b));return out
    endpoint_coeff=defaultdict(F);active=set();moving_gaps=0;max_active=0;trace=[]
    for j,(_,_,e) in enumerate(xrows):
        for sign,i in events[j]:
            if sign>0:
                if i in active:raise AssertionError('double insertion')
                active.add(i)
            else:
                if i not in active:raise AssertionError('unmatched deletion')
                active.remove(i)
        if j+1==len(xrows):break
        dv=xrows[j+1][2][1]-e[1]
        if not dv:continue
        moving_gaps+=1;max_active=max(max_active,len(active))
        unions=merge(ybounds[i] for i in active)
        for a,b in unions:
            endpoint_coeff[yrows[b][2]]+=dv;endpoint_coeff[yrows[a][2]]-=dv
        trace.append([j,dv,len(active),unions])
    if active:raise AssertionError('unclosed sweep')
    endpoint_coeff={e:c for e,c in endpoint_coeff.items() if c}
    constant=F(0);terms=defaultdict(F)
    for (a,v,s,q,b),c in endpoint_coeff.items():
        constant+=c*a
        if s:terms[q,b]+=c*s
    terms={key:c for key,c in terms.items() if c}
    beta=root_sum_eval(constant,terms,bits)
    betat=root_sum_eval(constant,terms,bits,derivative=True)
    beta_cell=root_sum_eval(constant,terms,bits,-tr,tr)
    betat_cell=root_sum_eval(constant,terms,bits,-tr,tr,derivative=True)
    M,Mp=map(F,[slab['uniform_mass']['center_total_mass'],slab['uniform_mass']['total_derivative']]);lam=F(slab['lambda'])
    if M<=abs(Mp)*tr:raise ValueError('mass denominator not uniformly positive')
    mixed=interval_scale(interval_sub(interval_scale(betat,M),interval_scale(beta,Mp)),lam/(M*M))
    mcell=(M-abs(Mp)*tr,M+abs(Mp)*tr)
    numerator=interval_sub(interval_mul(betat_cell,mcell),interval_scale(beta_cell,Mp))
    mixed_cell=interval_scale(interval_mul(numerator,(1/(mcell[1]*mcell[1]),1/(mcell[0]*mcell[0]))),lam)
    termlist=[{'q0':q,'q1':b,'coefficient':c} for (q,b),c in sorted(terms.items())]
    out={'status':'passed','scope':'exact joint raw-endpoint order and direct y-union root-sum derivatives; conditional on supplied bbox completeness/empty, source-catalogue and uniform-tail certificate',
      'source_candidate_code_sha256':slab['code_sha256'],'source_center_data_canonical_sha256':slab['center_data_canonical_sha256'],
      'base_dz_radius':Z,'base_dt_radius':T,'endpoint_rho':rho,'dz_radius':zr,'dt_radius':tr,
      'z0':slab['z0'],'t0':slab['t0'],'bits':bits,'endpoint_limiter':limiter,
      'x_endpoint_count':len(xrows),'y_endpoint_count':len(yrows),'consecutive_order_checks':checked,
      'all_center_equalities_persistent':True,'all_candidate_clipped_rectangles_positive':True,
      'box_count':len(boxes),'moving_x_gaps':moving_gaps,'max_active_at_moving_gap':max_active,
      'sweep_method':'direct sorted interval rank merge; no production UnionTree or endpoint sweep',
      'x_endpoint_order_sha256':digest([row[2] for row in xrows]),'y_endpoint_order_sha256':digest([row[2] for row in yrows]),
      'moving_gap_union_trace_sha256':digest(trace),'nonzero_y_endpoint_coefficients':[[e,c] for e,c in sorted(endpoint_coeff.items())],
      'beta_formula':'constant + sum coefficient*sqrt(q0+q1*(t-t0)); independent of z on certified rectangle',
      'beta_constant':constant,'beta_root_terms':termlist,'beta_center_bounds':beta,'beta_t_center_bounds':betat,
      'beta_cell_bounds':beta_cell,'beta_t_cell_bounds':betat_cell,
      'mass':M,'mass_t_derivative':Mp,'lambda':lam,
      'mixed_J_z_t_formula':'lambda*(M(t)*beta_t(t)-M_t*beta(t))/M(t)^2',
      'mixed_J_z_t_center_bounds':mixed,'mixed_J_z_t_cell_bounds':mixed_cell,
      'evidence_boundary':'No bbox enumeration or rejected-box loop rerun. Endpoint radius uses every raw active x/y endpoint and clipping boundary, so it is sufficient, not a maximal effective-union cell.',
      'python':sys.version,'platform':platform.platform(),'cpu_seconds':time.process_time()-start}
    with localcontext() as context:
        context.prec=30
        dec=lambda q:str(Decimal(q.numerator)/Decimal(q.denominator))
        out['decimals_display_only']={k:[dec(q) for q in out[k]] for k in ('beta_center_bounds','beta_t_center_bounds','mixed_J_z_t_center_bounds','mixed_J_z_t_cell_bounds')}
        out['decimals_display_only'].update({k:dec(out[k]) for k in ('dz_radius','dt_radius')})
    return out

def selftest():
    # Raw endpoint equality identity must retain the derivative of a square q0.
    a=expr([0,0,1,4,1]);b=expr([2,0,0,0,0])
    assert a!=b and center_key(a)==center_key(b)
    try:ordered({a,b},80)
    except ValueError:pass
    else:raise AssertionError('nonpersistent square center tie accepted')
    assert expr([0,0,1,F(4,9),0])==expr([F(2,3),0,0,0,0])
    for q in [F(1,9),F(2),F(4),F(1,10**100)]:
        lo,hi=root(q,160);assert lo*lo<=q<=hi*hi
    # One rectangle: x=[0,z+2], y=[0,sqrt(4+t)], A=(z+2)sqrt(4+t).
    e=lambda a,v=0,s=0,q=0,b=0:[str(a),str(v),s,str(q),str(b)]
    slab={'schema':'position-mass-slab-v1','status':'complete_candidate_and_empty_certificate',
      'dz_radius':'1/10','dt_radius':'1/10','domain':['-1','4','-1','4'],
      'boxes':[[e(0),e(2,1),e(0),e(0,0,1,4,1)]],
      'box_records':[{'bbox':[0,0,0,0]}],'counts':{'candidate_boxes':1},'z0':'0','t0':'0',
      'uniform_mass':{'center_total_mass':'1','total_derivative':'0'},'lambda':'1',
      'code_sha256':'toy','center_data_canonical_sha256':'toy'}
    out=certify(slab,160)
    assert out['beta_constant']==0 and out['beta_root_terms']==[{'q0':F(4),'q1':F(1),'coefficient':F(1)}]
    assert out['beta_center_bounds']==(F(2),F(2)) and out['beta_t_center_bounds']==(F(1,4),F(1,4))
    assert out['mixed_J_z_t_center_bounds']==(F(1,4),F(1,4))
    return {'status':'passed','tests':['square-center moving radical retained','nonpersistent equality rejected','constant rational-square normalization','exact dyadic square-root inclusions','one-rectangle beta and beta_t equal analytic 2 and 1/4'], 'toy_certificate':out}

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--slab',type=Path);ap.add_argument('--output',required=True,type=Path);ap.add_argument('--bits',type=int,default=160);ap.add_argument('--self-test',action='store_true');args=ap.parse_args()
    if args.self_test:out=selftest()
    elif args.slab:
        raw=args.slab.read_bytes();payload=gzip.decompress(raw) if args.slab.suffix=='.gz' else raw
        out=certify(json.loads(payload),args.bits);out['slab_file_sha256']=hashlib.sha256(raw).hexdigest();out['slab_payload_sha256']=hashlib.sha256(payload).hexdigest();out['slab_path_at_run']=str(args.slab.resolve())
    else:ap.error('--slab or --self-test required')
    out['script_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest();out['actual_command']=[sys.executable,*sys.argv]
    if not os.path.ismount('/Volumes/data'):raise RuntimeError('/Volumes/data is not mounted')
    args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(serial(out),indent=2)+'\n')
    print(json.dumps(serial({k:v for k,v in out.items() if k in ('status','scope','box_count','endpoint_rho','dz_radius','dt_radius','beta_center_bounds','beta_t_center_bounds','mixed_J_z_t_center_bounds','decimals_display_only','cpu_seconds','script_sha256')}),indent=2))
if __name__=='__main__':main()
