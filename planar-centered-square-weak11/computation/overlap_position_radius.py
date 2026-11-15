#!/usr/bin/env python3
"""R47 position radius from y-overlap-coupled x-endpoint comparisons.

Only velocities 0/1 are supported. For each outward-dyadic clipped y support,
activate both raw x endpoints, with the two F_x boundaries active throughout
F_y. Every positive-height y strip keeps its full active x endpoint order.
A count-aware segment tree finds the least directed positive adjacent gap
whose velocities differ. All coactive differing-velocity pairs are protected:
the path between any such pair contains an adjacent differing-velocity pair,
and its positive gap is no larger than the original pair's full gap.

The input complete enumeration and omitted-box safety certificates are an
external premise. This is a radius attachment, not a new complete enumeration
or an independent slope computation. The old global radius is diagnostic only.
"""
import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import sys
import time
from fractions import Fraction as F
from functools import cmp_to_key


def _serial(x):
    if isinstance(x,F):return str(x)
    if isinstance(x,dict):return {k:_serial(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):return [_serial(v) for v in x]
    return x


def _digest(x):
    return hashlib.sha256(json.dumps(_serial(x),sort_keys=True,separators=(',',':')).encode()).hexdigest()


def _load(path):
    path=Path(path).resolve();spec=importlib.util.spec_from_file_location('overlap_frozen_sweep',path)
    if spec is None or spec.loader is None:raise ValueError('cannot load frozen sweep')
    module=importlib.util.module_from_spec(spec);sys.path.insert(0,str(path.parent))
    try:spec.loader.exec_module(module)
    finally:sys.path.pop(0)
    return module


class CoupledGapTree:
    """Counted active endpoint classes, first/last and minimum mixed seam gap."""
    def __init__(self,order,bounds):
        self.order=order;self.bounds=bounds;self.n=len(order)
        self.count=[0]*self.n;size=4*self.n+4
        self.first=[None]*size;self.last=[None]*size;self.best=[None]*size
        self.mixed_seam_checks=0
    def update(self,index,delta):
        self.count[index]+=delta
        if self.count[index]<0:raise ValueError('negative canonical endpoint multiplicity')
        def change(node,lo,hi):
            if hi-lo==1:
                self.first[node]=self.last[node]=index if self.count[index] else None
                self.best[node]=None;return
            mid=(lo+hi)//2
            if index<mid:change(node*2,lo,mid)
            else:change(node*2+1,mid,hi)
            a,b=node*2,node*2+1
            self.first[node]=self.first[a] if self.first[a] is not None else self.first[b]
            self.last[node]=self.last[b] if self.last[b] is not None else self.last[a]
            options=[x for x in (self.best[a],self.best[b]) if x is not None]
            i,j=self.last[a],self.first[b]
            if i is not None and j is not None and self.order[i][1]!=self.order[j][1]:
                self.mixed_seam_checks+=1
                gap=self.bounds[j][0]-self.bounds[i][1]
                if gap<=0:raise ValueError('coactive mixed endpoint gap unresolved/nonpositive')
                options.append((gap,i,j))
            self.best[node]=min(options) if options else None
        change(1,0,self.n)
    def minimum(self):return self.best[1]


def overlap_radius(boxes,domain,H,delta_empty,*,sweep_module,bits=80,source_metadata=None):
    """Return new and old rational radii for a complete candidate raw box family.

    Include ALL original active boxes, even boxes empty after a recenter.
    H/delta_empty refer to the remaining original source/omitted-box certificates,
    and must not already be capped by the old global raw endpoint radius.
    """
    start=time.perf_counter();H,delta_empty=F(H),F(delta_empty)
    if min(H,delta_empty)<=0 or bits<1:raise ValueError('positive input radii/bits required')
    domain=tuple(map(F,domain));fl,fr,fb,ft=domain
    if fl>=fr or fb>=ft:raise ValueError('positive domain required')
    path=Path(sweep_module).resolve();sw=_load(path)
    left,right=sw.expr(fl),sw.expr(fr)
    family=[];allx={left,right};refs={left:[('domain','left')],right:[('domain','right')]}
    supports=[];permanent_y_empty=[]
    for index,box in enumerate(boxes):
        if len(box)!=4:raise ValueError('four raw endpoints per box required')
        l,r,y0,y1=(sw.expr(*e) for e in box)
        if l[1] not in (0,1) or r[1] not in (0,1):raise ValueError('only x velocities 0/1 supported')
        if y0[1] or y1[1]:raise ValueError('fixed y endpoints required')
        family.append((l,r,y0,y1));allx.update((l,r))
        refs.setdefault(l,[]).append((index,'left'));refs.setdefault(r,[]).append((index,'right'))
        a,_=sw.interval(y0,bits);_,b=sw.interval(y1,bits)
        lo,hi=max(fb,a),min(ft,b)
        if lo<hi:supports.append((index,lo,hi))
        else:permanent_y_empty.append({'box_index':index,'outer_clipped_support':[lo,hi],
                                      'reason':'outer dyadic support has nonpositive length; true y support permanently empty'})
    # This initial implementation fails closed at any nonpersistent central x
    # equality, including equalities of disjoint y supports. Crossing them away
    # from the center is permitted when no coupled comparison requires them.
    order=sorted(allx,key=cmp_to_key(lambda a,b:sw.compare(a,b,bits)))
    rank={e:i for i,e in enumerate(order)};bounds=[sw.interval(e,bits) for e in order]
    global_gap=None;global_pair=None
    for i,(a,b) in enumerate(zip(order,order[1:])):
        gap=bounds[i+1][0]-bounds[i][1]
        if gap<=0:raise ValueError('global center ordering not certified at requested bits')
        if a[1]!=b[1] and (global_gap is None or gap<global_gap):global_gap,global_pair=gap,(i,i+1)
    events={fb:[],ft:[]};support_by_index={}
    for index,lo,hi in supports:
        events.setdefault(lo,[]).append((index,1));events.setdefault(hi,[]).append((index,-1))
        support_by_index[index]=(lo,hi)
    tree=CoupledGapTree(order,bounds)
    tree.update(rank[left],1);tree.update(rank[right],1)
    active=set();best=None;limiter=None;positive_strips=0;event_updates=0
    levels=sorted(events)
    for k,y in enumerate(levels):
        # No zero-height transient strip contributes: all same-y events are
        # applied together before inspecting the next positive-height strip.
        for index,delta in events[y]:
            if delta==1:
                if index in active:raise ValueError('duplicate y enter')
                active.add(index)
            else:
                if index not in active:raise ValueError('y exit without enter')
                active.remove(index)
            l,r=family[index][:2];tree.update(rank[l],delta);tree.update(rank[r],delta)
            event_updates+=1
        if k+1==len(levels):continue
        upper=levels[k+1]
        if y>=upper:raise AssertionError('nonpositive y strip')
        positive_strips+=1;candidate=tree.minimum()
        if candidate is not None and (best is None or candidate<best):
            best=candidate;gap,i,j=candidate
            def owner(endpoint):
                for index,side in refs[endpoint]:
                    if index=='domain' or index in active:return [index,side]
                raise AssertionError('limiting endpoint has no active owner')
            limiter={'endpoint_ranks':[i,j],'endpoints':[order[i],order[j]],
                     'gap_lower':gap,'relative_speed':1,'shared_outer_y_strip':[y,upper],
                     'active_box_count_in_strip':len(active),'endpoint_owners':[owner(order[i]),owner(order[j])]}
    if event_updates!=2*len(supports):raise AssertionError('incomplete y support event ledger')
    if active:raise AssertionError('unclosed y support events')
    expected=[0]*len(order);expected[rank[left]]+=1;expected[rank[right]]+=1
    if tree.count!=expected:raise AssertionError('canonical multiplicities not restored after sweep')
    coupled_gap=best[0] if best is not None else None
    coupled_radius=coupled_gap/2 if coupled_gap is not None else None
    old_radius=global_gap/2 if global_gap is not None else None
    constraints=[('remaining_source_H',H),('remaining_original_empty_radius',delta_empty)]
    if coupled_radius is not None:constraints.append(('overlap_coupled_endpoint_order',coupled_radius))
    reason,final=min(constraints,key=lambda z:z[1])
    old_final=min([H,delta_empty]+([old_radius] if old_radius is not None else []))
    if final<old_final:raise AssertionError('overlap radius unexpectedly below global directed radius')
    coverage={
        'candidate_boxes':len(family),'positive_outer_y_supports':len(supports),
        'permanently_y_empty':len(permanent_y_empty),'all_raw_x_endpoint_classes':len(order),
        'y_levels':len(levels),'positive_height_y_strips':positive_strips,'box_enter_exit_events':event_updates,
        'canonical_endpoint_count_updates':2*event_updates+2,'segment_tree_mixed_seam_checks':tree.mixed_seam_checks,
        'velocity_domain_validated':[0,1],'same_y_events_grouped':True,'all_candidate_raw_endpoints_retained':True,
        'Fx_boundaries_active_throughout_Fy':True,'canonical_multiplicity_restored':True,
        'proof':'Every pair with true positive y overlap is coactive on a positive outer-support strip. On that strip any mixed-velocity pair path contains an adjacent mixed-velocity pair with no larger true gap. All such adjacent pairs are maintained by first/last/min segment-tree merges. Their directed gap lower bounds divided by two preserve every coactive order for velocities 0/1. Each box has its own two raw endpoints and both F_x endpoints coactive wherever its true y support is nonempty, preserving self length and clipping. Fixed strip interval orders make the full union area affine; original omitted boxes remain empty by delta_empty.',
        'scope':'finite complete candidate family supplied by independently accepted input; no new general C2 bound'}
    out={'schema':'overlap-position-radius-v1','status':'certified_radius_attachment','bits':bits,
         'mathematical_contract_sha256':'3e391beace80ec39a668eaeb81de218d8ab2d982da39d919f8f3ead749d5a3b7',
         'external_premise':'The complete original candidate family, source-stable H, omitted-box delta_empty, and position-independent tail were independently accepted by caller. No old global endpoint radius is used as a final cap.',
         'domain':list(domain),'H':H,'delta_empty':delta_empty,'old_global_raw_radius':old_radius,
         'old_global_final_delta':old_final,'old_global_limiting_pair':global_pair,
         'overlap_gap_lower':coupled_gap,'overlap_radius':coupled_radius,'final_delta':final,
         'final_delta_limiter':reason,'radius_constraints':constraints,'limiting_pair':limiter,
         'improvement_factor':final/old_final,'coverage':coverage,'endpoint_order':order,
         'outer_y_support_catalog':supports,'permanently_y_empty_boxes':permanent_y_empty,
         'source_metadata':source_metadata or {},'input_candidate_boxes_sha256':_digest(boxes),
         'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
         'sweep_module_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
         'interval_primitives_sha256':hashlib.sha256(path.with_name('r039_interval_primitives.py').read_bytes()).hexdigest(),
         'python':sys.version,'elapsed_seconds':time.perf_counter()-start}
    return _serial(out)


def _brute_radius(out):
    # Independent quadratic coactive-pair oracle for small tests. Unlike the
    # production tree this explicitly checks every pair on every y strip.
    order=[tuple((F(e[0]),F(e[1]),e[2],F(e[3]))) for e in out['endpoint_order']]
    # Selftests use rational endpoints, so these exact values need no root tool.
    assert all(e[2]==0 for e in order)
    domain=list(map(F,out['domain']));levels={domain[2],domain[3]}
    supports=out['outer_y_support_catalog']
    for _,a,b in supports:levels.update((F(a),F(b)))
    best=None;family=_BRUTE_TEST_FAMILY
    for y,z in zip(sorted(levels),sorted(levels)[1:]):
        endpoints={(domain[0],F(0)),(domain[1],F(0))}
        for index,a,b in supports:
            if F(a)<=y and z<=F(b):
                for e in family[index][:2]:endpoints.add((F(e[0]),F(e[1])))
        endpoints=sorted(endpoints)
        for i,(a,va) in enumerate(endpoints):
            for b,vb in endpoints[i+1:]:
                if va!=vb:
                    gap=b-a
                    if gap<=0:raise AssertionError('oracle central collision')
                    if best is None or gap<best:best=gap
    return None if best is None else best/2


def selftest(sweep_module):
    global _BRUTE_TEST_FAMILY
    def e(a,v=0):return [str(a),str(v),0,'0']
    def rect(a,b,c,d,vl=0,vr=0):return [e(a,vl),e(b,vr),e(c),e(d)]
    domain=['-10','10','-4','4'];cases=[]
    def check(name,family):
        global _BRUTE_TEST_FAMILY
        _BRUTE_TEST_FAMILY=family
        out=overlap_radius(family,domain,5,5,sweep_module=sweep_module)
        assert _brute_radius(out)==(None if out['overlap_radius'] is None else F(out['overlap_radius']))
        cases.append({'name':name,'old':out['old_global_raw_radius'],'new':out['overlap_radius'],
                      'final':out['final_delta'],'limiting_pair':out['limiting_pair']})
        return out
    disjoint=check('disjoint_support_x_crossing_ignored',[rect(0,2,0,1),rect('201/100','301/100',2,3,1,1)])
    assert F(disjoint['overlap_radius'])>F(disjoint['old_global_raw_radius'])
    overlap=check('positive_y_overlap_restricts',[rect(0,2,0,2),rect('201/100','301/100',1,3,1,1)])
    assert overlap['overlap_radius']==overlap['old_global_raw_radius']=='1/200'
    touch=check('same_y_exit_enter_only_touch',[rect(0,2,0,1),rect('201/100','301/100',1,2,1,1)])
    assert F(touch['overlap_radius'])>F(touch['old_global_raw_radius'])
    counted=check('persistent_endpoint_multiplicity',[rect(0,2,0,1),rect(0,2,0,2),rect('201/100','301/100','3/2','5/2',1,1)])
    assert counted['overlap_radius']=='1/200'
    clipped=check('Fx_clip_branch_restricts',[rect('999/100','1099/100',0,1,1,1)])
    assert clipped['overlap_radius']=='1/200'
    # The true y supports only touch at sqrt(2). Outward dyadic rounding
    # creates a false overlap, which must conservatively retain the x event.
    radical_touch=check('outward_radical_y_false_overlap',[
        [e(0),e(1),e(0),['0','0',1,'2']],
        [e('1001/1000',1),e('2001/1000',1),['0','0',1,'2'],e(3)]])
    assert radical_touch['overlap_radius']==radical_touch['old_global_raw_radius']=='1/2000'
    # Persistent clipping equality shares a canonical class with F_x.
    persistent_clip=check('persistent_Fx_endpoint_equality',[rect(10,11,0,1)])
    assert persistent_clip['overlap_radius'] is None
    radical_x=overlap_radius([
        [e(0),['0','0',1,'2'],e(0),e(1)],
        [['1/1000','1',1,'2'],e(3,1),e(2),e(3)]],domain,5,5,sweep_module=sweep_module)
    assert F(radical_x['overlap_radius'])>F(radical_x['old_global_raw_radius'])
    import random
    rng=random.Random(47001)
    for case in range(32):
        family=[]
        for i in range(5):
            a=F(rng.randrange(-5,5))+F(i+1,1000);b=a+rng.randrange(1,4)
            c=rng.randrange(-3,3);d=c+rng.randrange(1,4)
            family.append(rect(a,b,c,d,rng.randrange(2),rng.randrange(2)))
        check('random_'+str(case),family)
    try:overlap_radius([rect(0,1,0,1,0,2)],domain,1,1,sweep_module=sweep_module)
    except ValueError:pass
    else:raise AssertionError('nonbinary velocity accepted')
    try:overlap_radius([rect(0,1,0,1),rect(1,2,2,3,1,1)],domain,1,1,sweep_module=sweep_module)
    except ValueError:pass
    else:raise AssertionError('nonpersistent central equality accepted')
    return {'status':'passed','seed':47001,'quadratic_oracle_cases':len(cases),'cases':cases,
            'radical_x_disjoint_improvement':{'old':radical_x['old_global_raw_radius'],'new':radical_x['overlap_radius']},
            'nonbinary_velocity_rejected':True,'nonpersistent_central_equality_rejected':True,
            'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('input',nargs='?')
    ap.add_argument('--sweep-module',required=True);ap.add_argument('--bits',type=int,default=80)
    ap.add_argument('--accepted-input-sha256');ap.add_argument('--self-test',action='store_true');ap.add_argument('--output',required=True)
    args=ap.parse_args()
    if args.self_test:out=selftest(args.sweep_module)
    else:
        if not args.input:ap.error('input JSON required')
        path=Path(args.input);digest=hashlib.sha256(path.read_bytes()).hexdigest()
        if args.accepted_input_sha256 and digest!=args.accepted_input_sha256:raise ValueError('accepted input SHA256 mismatch')
        data=json.loads(path.read_text());schema=data.get('schema')
        if schema=='position-slab-enumerator-v1':
            boxes=data['boxes'];H=data['H'];empty=data['delta_empty']
        elif schema=='recenter-position-slab-v1':
            boxes=data['candidate_boxes'];H=data['remaining_source_H'];empty=data['remaining_original_empty_radius']
        else:raise ValueError('requires original complete slab or full-candidate recenter attachment')
        metadata={k:data[k] for k in ('schema','s0','code_sha256','source_catalog_sha256','input_source_catalog_sha256') if k in data}
        out=overlap_radius(boxes,data['domain'],H,empty,sweep_module=args.sweep_module,bits=args.bits,source_metadata=metadata)
        out['input_file_sha256']=digest;out['command']=[sys.executable,*sys.argv]
    if not os.path.ismount('/Volumes/data'):raise RuntimeError('/Volumes/data is not mounted')
    target=Path(args.output);target.parent.mkdir(parents=True,exist_ok=True);target.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k in ('schema','status','old_global_raw_radius','overlap_radius','final_delta','improvement_factor','coverage','code_sha256','quadratic_oracle_cases')},indent=2))


if __name__=='__main__':main()
