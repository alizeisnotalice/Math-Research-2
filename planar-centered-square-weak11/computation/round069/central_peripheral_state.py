#!/usr/bin/env python3
r"""Exact W/Wprime, internal-halo and pure-mixed K/periphery diagnostics.

Frozen root/rectangle libraries are reused. Only a rational integer dilation
clears coordinate denominators; the source domain is never normalized by
sqrt(L) and no algebraic source coordinates are introduced.
"""
import argparse
from fractions import Fraction as F
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import sys
import time


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def load(name,p):
    s=importlib.util.spec_from_file_location(name,p);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--archive-root',type=Path,required=True)
    ap.add_argument('--spec',type=Path,required=True);ap.add_argument('--job',required=True)
    ap.add_argument('--out',type=Path,required=True);args=ap.parse_args();assert os.path.ismount('/Volumes/data')
    sp=args.spec if args.spec.is_absolute() else args.archive_root/args.spec
    spec=json.loads(sp.read_text());job=next(j for j in spec['jobs'] if j['job_id']==args.job)
    ip=args.archive_root/job['input'];assert sha(ip)==job['input_sha256']
    gp=args.archive_root/spec['geometry_source'];ep=args.archive_root/spec['rectclip_source']
    assert sha(gp)==spec['geometry_sha256'];assert sha(ep)==spec['rectclip_sha256']
    mod=load('frozen_r064_geometry',gp);er=load('frozen_r066_rectclip',ep)
    started=time.perf_counter();data=json.loads(ip.read_text());points=[tuple(map(F,p)) for p in data['points']]
    mass=list(map(F,data['initial_masses']));n=len(points);S=data['S_mask'];assert 0<S<(1<<n)-1
    assert min(mass)>0 and len(points)==len(set(points))==len(mass)
    c=math.lcm(*(v.denominator for p in points for v in p));bits=spec['bits'];D0=1<<bits
    D=c*D0;den=D*D;ps=[(x*c,y*c) for x,y in points];ms=[m*c*c for m in mass]
    assert all(v.denominator==1 for p in ps for v in p)
    boxes=mod.closures(ps);eligible=[b for b in boxes if sum(ms[i] for i in b['ids'])>b['span2']]
    maxscaled=max(sum(ms[i] for i in b['ids']) for b in eligible)
    maxima=sorted(b['mask'] for b in eligible if sum(ms[i] for i in b['ids'])==maxscaled)
    L=sum(mass[i] for i in range(n) if S>>i&1);Ls=L*c*c
    assert maxscaled==Ls and S in maxima and maxima==data['expected_maximum_masks']
    assert L==F(data['expected_L'])
    full_i,full_o=mod.exact_rectangles(ps,ms,boxes,bits)
    f_i={r[4]:[r] for r in full_i};f_o={r[4]:[r] for r in full_o}
    wi=[];wo=[];pi=[];po=[];labels=[]
    for B in eligible:
        t=B['mask']
        if not t&S:continue
        outside=[i for i in B['ids'] if not S>>i&1]
        bm=sum((ms[i] for i in outside),F(0));a=sum(ms[i] for i in B['ids'])-bm
        ki,ko=er.receiver(B['bbox'],bm,bits,t)
        if outside:
            tight=(min(ps[i][0] for i in outside),max(ps[i][0] for i in outside),
                   min(ps[i][1] for i in outside),max(ps[i][1] for i in outside))
            ti,to=er.receiver(tight,bm,bits,t)
        else:ti,to=[],[]
        ii=f_i.get(t,[]);oo=f_o.get(t,[])
        wi += [q for r in ii for q in er.subtract_rect(r,ko)]
        wo += [q for r in oo for q in er.subtract_rect(r,ki)]
        pi += [q for r in ii for q in er.subtract_rect(r,to)]
        po += [q for r in oo for q in er.subtract_rect(r,ti)]
        labels.append({'mask':t,'rank':t.bit_count(),'mass':str((a+bm)/(c*c)),
                       'a':str(a/(c*c)),'b':str(bm/(c*c)),
                       'full_bbox':list(map(lambda v:str(v/c),B['bbox'])),
                       'full_span_squared':str(B['span2']/(c*c))})
    ids=[i for i in range(n) if S>>i&1];keep=[i for i in range(n) if not S>>i&1]
    ins=[ps[i] for i in ids];ims=[ms[i] for i in ids]
    rem=[ps[i] for i in keep];rms=[ms[i] for i in keep]
    ib=mod.closures(ins);rb=mod.closures(rem)
    in_i,in_o=mod.exact_rectangles(ins,ims,ib,bits);rem_i,rem_o=mod.exact_rectangles(rem,rms,rb,bits)
    a=min(ps[i][0] for i in ids);b=max(ps[i][0] for i in ids)
    y0=min(ps[i][1] for i in ids);y1=max(ps[i][1] for i in ids)
    ri,ro=er.radii(3*Ls,bits)
    def central(r):return (int(b*D0)-r,int(a*D0)+r,int(y1*D0)-r,int(y0*D0)+r,0)
    Ki=central(ri);Ko=central(ro);assert Ki[0]<Ki[1] and Ki[2]<Ki[3]
    # Area memoization ignores source-mask annotations, which do not change
    # unions. Every rectangle is still processed with exact integer endpoints.
    cache={}
    def canonical(rs):return tuple(sorted({tuple(r[:4]) for r in rs}))
    def DA(A,B):
        x,y=canonical(A),canonical(B);key=(x,y)
        if key not in cache:
            cache[key]=0 if not x or x==y else mod.difference_area([r+(0,) for r in x],[r+(0,) for r in y])
        return cache[key]
    def UA(A):return DA(A,[])
    def clip(A,K):
        out=[]
        for r in A:
            q=(max(r[0],K[0]),min(r[1],K[1]),max(r[2],K[2]),min(r[3],K[3]),0)
            if q[0]<q[1] and q[2]<q[3]:out.append(q)
        return out
    def outside(A,K):return [q for r in A for q in er.subtract_rect(r,[K])]
    def enc(lo,hi,signed=False):
        assert lo<=hi and (signed or lo>=0),(lo,hi)
        return {'lower':str(F(lo,den)),'upper':str(F(hi,den)),
                'normalized_lower':str(F(lo,den)/L),'normalized_upper':str(F(hi,den)/L)}
    regions={};numbers={}
    def region(name,kind,Ai,Ao,Bi=None,Bo=None,Ci=None,Co=None):
        values={}
        for part,li,uo in [('global',Ai,Ao),('central',clip(Ai,Ki),clip(Ao,Ko)),
                            ('peripheral',outside(Ai,Ko),outside(Ao,Ki))]:
            if kind=='union':lo,hi=UA(li),UA(uo)
            elif kind=='difference':lo,hi=DA(li,Bo),DA(uo,Bi)
            elif kind=='intersection':lo,hi=UA(li)-DA(li,Bi),UA(uo)-DA(uo,Bo)
            elif kind=='intersection_minus':
                lo=DA(li,Co)-DA(li,Bi+Co);hi=DA(uo,Ci)-DA(uo,Bo+Ci)
            else:raise ValueError(kind)
            assert 0<=lo<=hi;values[part]=(lo,hi)
        g=values['global'];k=values['central'];o=values['peripheral']
        residual=(g[0]-k[1]-o[1],g[1]-k[0]-o[0]);assert residual[0]<=0<=residual[1]
        numbers[name]=values;regions[name]={part:enc(*v) for part,v in values.items()}
        regions[name]['partition_identity_residual']=enc(*residual,signed=True)
    region('Hfull','union',full_i,full_o)
    region('R','union',rem_i,rem_o)
    region('I','union',in_i,in_o)
    region('W','union',wi,wo)
    region('Wprime','union',pi,po)
    region('Delta','difference',full_i,full_o,rem_i,rem_o)
    region('P_pure_mixed','difference',full_i,full_o,in_i+rem_i,in_o+rem_o)
    region('internal_old_credit','intersection',in_i,in_o,rem_i,rem_o)
    region('W_old_surplus','intersection',wi,wo,rem_i,rem_o)
    region('Wprime_old_surplus','intersection',pi,po,rem_i,rem_o)
    region('Wprime_surplus_outside_I','intersection_minus',pi,po,rem_i,rem_o,in_i,in_o)
    checks={'I_outer_minus_Hfull_outer':DA(in_o,full_o),'R_outer_minus_Hfull_outer':DA(rem_o,full_o),
            'Wprime_inner_minus_W_inner':DA(pi,wi),'Wprime_outer_minus_W_outer':DA(po,wo),
            'W_outer_minus_Hfull_outer':DA(wo,full_o),
            'Hfull_inner_minus_Hrem_outer_union_Wprime_inner':DA(full_i,rem_o+pi)}
    assert all(v==0 for v in checks.values()),checks
    identities={};net={}
    for part in ['global','central','peripheral']:
        dl,du=numbers['Delta'][part];il,iu=numbers['I'][part]
        cl,cu=numbers['internal_old_credit'][part];pl,pu=numbers['P_pure_mixed'][part]
        residual=(dl-(iu-cl+pu),du-(il-cu+pl));assert residual[0]<=0<=residual[1]
        identities[f'Delta_minus_rhs_{part}']=enc(*residual,signed=True)
        nlo,nhi=pl-cu,pu-cl
        net[part]={'P_minus_credit':enc(nlo,nhi,signed=True),
                   'P_at_most_L_over_5_certified':F(pu,den)<=L/5,
                   'P_strictly_above_L_over_5_certified':F(pl,den)>L/5,
                   'net_at_most_L_over_5_certified':F(nhi,den)<=L/5,
                   'net_strictly_above_L_over_5_certified':F(nlo,den)>L/5}
        for state,surplus in [('W','W_old_surplus'),('Wprime','Wprime_old_surplus')]:
            wl,wu=numbers[state][part];sl,su=numbers[surplus][part]
            r=(wl-du-su,wu-dl-sl);assert r[0]<=0<=r[1]
            identities[f'{state}_minus_rhs_{part}']=enc(*r,signed=True)
    rcentral=numbers['R']['central'];assert rcentral[0]==0
    lo=max(numbers[name]['central'][0] for name in ['Hfull','Delta','W','Wprime'])
    hi=min(numbers[name]['central'][1] for name in ['Hfull','Delta','W','Wprime']);assert lo<=hi
    component_masks=[1<<i for i in range(n)]
    for B in eligible:
        mask=B['mask'];joined=mask
        for cc in component_masks:
            if cc&mask:joined |= cc
        component_masks=[cc for cc in component_masks if not cc&mask]+[joined]
    total=sum(mass);span=max(max(p[t] for p in points)-min(p[t] for p in points) for t in (0,1))
    output={'schema':'r069-central-peripheral-state-v1','status':'PASS','id':data['id'],'source_input_id':data['source_input_id'],
      'python':sys.version,'source_sha256':sha(Path(__file__)),'spec_sha256':sha(sp),'input_sha256':sha(ip),
      'geometry_source_sha256':sha(gp),'rectclip_source_sha256':sha(ep),'source_count':n,'S_mask':S,'S_rank':S.bit_count(),'L':str(L),
      'all_maximum_masks':maxima,'complete_closed_bbox_count':len(boxes),'strict_eligible_count':len(eligible),
      'eligible_intersect_S_count':len(labels),'internal_closed_bbox_count':len(ib),'residual_closed_bbox_count':len(rb),
      'eligible_component_masks':component_masks,'whole_nonterminal_margin':str(span*span-F(35,33)*total),
      'S_bbox':[str(a/c),str(b/c),str(y0/c),str(y1/c)],
      'S_width_squared_over_L':[str(((b-a)/c)**2/L),str(((y1-y0)/c)**2/L)],
      'computational_integer_dilation':c,'bits_beyond_rational_source_grid':bits,'physical_grid_step':str(F(1,D)),
      'K':{'formula':'[xmax-rho,xmin+rho]x[ymax-rho,ymin+rho],rho=sqrt(3L)/2',
           'inner':list(map(lambda v:str(F(v,D)),Ki[:4])),'outer':list(map(lambda v:str(F(v,D)),Ko[:4])),
           'rho_lower':str(F(ri,D)),'rho_upper':str(F(ro,D)),'rho_squared':'3L/4'},
      'regions':regions,'area_identities':identities,
      'identity_formulas':{'Delta':'I-internal_old_credit+P_pure_mixed',
                           'W':'Delta+W_old_surplus','Wprime':'Delta+Wprime_old_surplus'},
      'mixed_mass_gap_diagnostics':net,
      'central_residual_halo_geometrically_certified_empty':rcentral[1]==0,
      'central_Hfull_Delta_W_Wprime_common_interval':enc(lo,hi),
      'exact_zero_area_inclusion_checks':{k:str(v) for k,v in checks.items()},
      'all_intersecting_labels':labels,'area_cache_entries':len(cache),'wall_seconds':time.perf_counter()-started,
      'scope':'Exactly one fixed maximum state. W/Wprime remain full unions without Hremaining discount; I/R/P/central partitions are separate diagnostics. No algebraic-coordinate normalization, theta scan, new random input, LP or universal .2 claim.'}
    args.out.parent.mkdir(parents=True,exist_ok=True);args.out.write_text(json.dumps(output,indent=2)+'\n')
    print(json.dumps({'status':'PASS','id':data['id'],'Wprime_upper':regions['Wprime']['global']['normalized_upper'],
                     'P_upper':regions['P_pure_mixed']['global']['normalized_upper'],
                     'net_upper':net['global']['P_minus_credit']['normalized_upper'],
                     'R_intersect_K_zero':rcentral[1]==0,'seconds':output['wall_seconds']}))


if __name__=='__main__':main()
