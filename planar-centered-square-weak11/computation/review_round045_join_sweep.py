#!/usr/bin/env python3
"""R45 only: reuse frozen R39 independent kernels; replay 11 new active catalogues.
No production sweep, candidate-bbox enumeration, or old tiny tests are run.
"""
import argparse,ast,gzip,hashlib,json,platform,sys,time
from decimal import Decimal,localcontext
from fractions import Fraction as F
from functools import lru_cache
from math import isqrt,lcm
from pathlib import Path

parser=argparse.ArgumentParser(description=__doc__)
for name in ['index','operations','reference-audit','prior-check','archive-repository','output']:
    parser.add_argument('--'+name,required=True,type=Path)
args=parser.parse_args()
for k,v in vars(args).items():setattr(args,k,v.resolve())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def json_conv(x):
    if isinstance(x,F):return str(x)
    if isinstance(x,dict):return {str(k):json_conv(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):return [json_conv(v) for v in x]
    return x
assert sha(args.index)=='855c330b1d36dc212956eefcfba840a831aeabe17d7f8ba6dfa94b5466986a88'
assert sha(args.reference_audit)=='46ac35f415025e03af461b9ea26a71f1c2a71c04968389ef1dbc66a6f4057bd2'
assert sha(args.prior_check)=='2667aa8d1803bfe51ff11bfaddd9a735405da51ac2a1c1a508b6be511ecd7219'
prior=json.loads(args.prior_check.read_text());P=F(prior['real']['normalization']['base_Px'])
source_sha={name:sha(args.operations/name) for name in prior['source_sha256']}
assert source_sha==prior['source_sha256']
source=args.reference_audit.read_text();tree=ast.parse(source)
kernels=[n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name in ['root','canon','enclosure','merge']]
assert len(kernels)==4
namespace={'F':F,'isqrt':isqrt,'lru_cache':lru_cache}
exec(compile(ast.Module(body=kernels,type_ignores=[]),str(args.reference_audit),'exec'),namespace)
canon,enclosure,merge=(namespace[k] for k in ['canon','enclosure','merge'])
index=json.loads(args.index.read_text());assert len(index['runs'])==11
assert len({x['name'] for x in index['runs']})==11
start=time.process_time();summary=[]
def loadraw(folder,name):
    p=folder/name
    if p.exists():return p.read_bytes(),p
    q=folder/(name+'.gz')
    return gzip.decompress(q.read_bytes()),q

def replay(entry):
    st=time.process_time();folder=args.index.parent/entry['name'];hashes={};data={}
    for name,h in entry['hashes'].items():
        raw,p=loadraw(folder,name);assert hashlib.sha256(raw).hexdigest()==h;hashes[name]=h
        data[name]=json.loads(raw)
    slab,slope,manifest,velocities=(data[k] for k in ['slab.json','slope.json','run_manifest.json','velocities.json'])
    assert slope['slab_sha256']==hashes['slab.json'] and slope['slab_code_sha256']==slab['code_sha256']==manifest['enum_sha256']
    assert slope['script_sha256']==source_sha['r039_join_slope.py']==manifest['sweep_sha256']
    for name,h in slope['dependency_sha256'].items():assert h==source_sha[name]
    assert 'conditional' in slope['scope'] and slab['status']=='complete'
    assert slab['velocities']==velocities and manifest['velocity_sha256']==hashes['velocities.json']
    assert manifest['input_sha256']==slab['input_file_sha256']
    assert F(slab['s0'])==F(slope['s0'])==F(manifest['physical_s0'])==P*F(entry['normalized_s0'])
    assert F(manifest['normalized_s0'])==F(entry['normalized_s0'])
    assert F(manifest['physical_H'])==P*F(manifest['normalized_H'])==P*F(entry['normalized_requested_H'])==F(slab['requested_H'])
    assert slab['counts']==entry['counts'] and slab['tail']==slope['tail']==entry['tail']
    bits=slope['sweep']['bits'];assert bits==80
    xl,xr,yl,yr=map(F,slab['domain']);assert xr-xl==2*P
    left,right=canon([xl,0,0,0]),canon([xr,0,0,0])
    boxes=[tuple(canon(e) for e in box) for box in slab['boxes']]
    assert len(boxes)==slab['counts']['active_boxes']==len(slab['box_records'])==slope['sweep']['box_count']
    allx={left,right}
    for l,r,y0,y1 in boxes:
        assert y0[1]==y1[1]==0;allx.update((l,r))
    rows=sorted((enclosure(e,bits)[0],enclosure(e,bits)[1],e) for e in allx)
    assert len(rows)==slope['sweep']['endpoint_count']
    rank={row[2]:i for i,row in enumerate(rows)};endpoint_radius=None;limiter=None
    for a,b in zip(rows,rows[1:]):
        gap=b[0]-a[1];assert gap>0;dv=abs(b[2][1]-a[2][1])
        if dv:
            d=gap/(2*dv)
            if endpoint_radius is None or d<endpoint_radius:endpoint_radius,limiter=d,(a[2],b[2])
    assert endpoint_radius==F(slope['sweep']['endpoint_radius'])
    assert limiter==tuple(canon(e) for e in slope['sweep']['limiting_pair'])
    events=[[] for _ in rows];ybounds=[];S=1
    for i,(l,r,y0,y1) in enumerate(boxes):
        lo=max(rank[l],rank[left]);hi=min(rank[r],rank[right]);assert lo<hi
        a,b=enclosure(y0,bits);c,d=enclosure(y1,bits)
        inner=max(yl,b),min(yr,c);outer=max(yl,a),min(yr,d)
        assert outer[0]<=inner[0]<inner[1]<=outer[1]
        ybounds.append((inner,outer));events[lo].append((1,i));events[hi].append((-1,i))
        for z in (*inner,*outer):S=lcm(S,z.denominator)
    yints=[tuple(tuple(int(z*S) for z in iv) for iv in pair) for pair in ybounds]
    active=set();beta_lo=beta_hi=F(0);gaps=0;max_active=0
    for j,row in enumerate(rows):
        for sign,i in events[j]:
            if sign>0:assert i not in active;active.add(i)
            else:assert i in active;active.remove(i)
        if j==len(rows)-1:break
        dv=rows[j+1][2][1]-row[2][1]
        if dv:
            gaps+=1;max_active=max(max_active,len(active))
            lo=F(merge(yints[i][0] for i in active),S);hi=F(merge(yints[i][1] for i in active),S)
            assert lo<=hi
            beta_lo+=dv*(lo if dv>0 else hi);beta_hi+=dv*(hi if dv>0 else lo)
    assert not active
    assert beta_lo==F(slope['beta_lower'])==F(slope['sweep']['beta_lower'])==F(entry['beta_lower'])
    assert beta_hi==F(slope['beta_upper'])==F(slope['sweep']['beta_upper'])==F(entry['beta_upper'])
    assert beta_lo<=beta_hi
    sign='positive' if beta_lo>0 else 'negative' if beta_hi<0 else 'unresolved'
    assert sign!='unresolved'
    H,empty,delta=map(F,[slab['H'],slab['delta_empty'],slope['delta']])
    assert delta==min(H,empty,endpoint_radius)==F(entry['physical_delta'])>0
    lam,mass=map(F,[slab['tail']['lambda'],slab['tail']['mass']]);assert F(slab['tail']['bound'])<lam and mass>0
    result={'name':entry['name'],'hashes':hashes,'input_sha256':manifest['input_sha256'],'enum_sha256':slab['code_sha256'],
        'velocity_nonzero_labels':[i for i,v in enumerate(velocities) if v],
        's0_physical':F(slab['s0']),'s0_normalized':F(slab['s0'])/P,'H':H,'delta_empty':empty,'endpoint_radius':endpoint_radius,
        'delta_physical':delta,'delta_normalized':delta/P,'beta_lower':beta_lo,'beta_upper':beta_hi,'sign':sign,
        'dJ_dphysical_bounds':(lam*beta_lo/mass,lam*beta_hi/mass),
        'dJ_dnormalized_bounds':(P*lam*beta_lo/mass,P*lam*beta_hi/mass),
        'boxes':len(boxes),'raw_endpoint_classes':len(rows),'nonzero_velocity_gaps_replayed':gaps,'max_active_at_replayed_gap':max_active,
        'tail_bound':F(slab['tail']['bound']),'tail_radius':F(slab['tail']['radius']),
        'all_raw_gaps_and_clips_checked':True,'cpu_seconds':time.process_time()-st}
    print(json.dumps({'run':entry['name'],'sign':sign,'boxes':len(boxes),'gaps':gaps,'cpu_seconds':result['cpu_seconds']}),flush=True)
    return result
for entry in index['runs']:summary.append(replay(entry))
lookup={r['name']:r for r in summary};pairs=[]
for direction,left_name,right_name in [('q0','local-slope-q0','negative-slope-q0'),('q1','local-slope','negative-slope-q1'),('q2','local-slope-q2','negative-slope-q2'),('q3','local-slope-q3','negative-slope-q3'),('all','local-slope-all','negative-slope-all-next')]:
    left,right=lookup[left_name],lookup[right_name]
    assert left['sign']=='positive' and right['sign']=='negative'
    assert left['s0_normalized']<right['s0_normalized']
    assert left['hashes']['velocities.json']==right['hashes']['velocities.json'] and left['input_sha256']==right['input_sha256']
    pairs.append({'direction':direction,'left_run':left_name,'right_run':right_name,
                  'open_normalized_bracket':[left['s0_normalized'],right['s0_normalized']],
                  'interior_maximum_or_plateau_exists_conditional_on_catalogue_audit':True,
                  'uniqueness_or_multidimensional_optimality_claimed':False})
mid=lookup['negative-slope-all'];assert mid['sign']=='positive'
assert lookup['local-slope-all']['s0_normalized']<mid['s0_normalized']<lookup['negative-slope-all-next']['s0_normalized']
archive_manifest=args.archive_repository/'planar-centered-square-weak11/computation/results/round045/source-manifest.json'
archive=json.loads(archive_manifest.read_text());archive_checks=[]
for row in archive['files']:
    p=args.archive_repository/row['destination'];raw=p.read_bytes();assert len(raw)==row['archived_size'] and hashlib.sha256(raw).hexdigest()==row['archived_sha256']
    uncompressed=gzip.decompress(raw) if p.suffix=='.gz' else raw
    assert len(uncompressed)==row['raw_size'] and hashlib.sha256(uncompressed).hexdigest()==row['raw_sha256']
    archive_checks.append({'destination':row['destination'],'raw_sha256':row['raw_sha256'],'archived_sha256':row['archived_sha256']})
arch_index=args.archive_repository/'planar-centered-square-weak11/computation/results/round045/results_index.json'
assert sha(arch_index)==sha(args.index)
result={'status':'passed_sweep_replay_pending_catalogue_binding','scope':'11 new active catalogues only; reused R39 proof/source tests, no old tests or bbox enumeration rerun',
    'index_sha256':sha(args.index),'source_sha256':source_sha,'reused_independent_script_sha256':sha(args.reference_audit),
    'reused_prior_check_sha256':sha(args.prior_check),'base_horizontal_period':P,'runs':summary,'five_direction_pairs':pairs,
    'extra_all_midpoint_positive':{'name':mid['name'],'s0_normalized':mid['s0_normalized'],'beta_lower':mid['beta_lower'],'beta_upper':mid['beta_upper']},
    'archive_manifest_sha256':sha(archive_manifest),'archive_files_checked':archive_checks,
    'not_recomputed':['candidate-bbox loops','rejected-box classifications','periodic halo areas','frozen R39 tiny tests'],
    'total_main_cpu_seconds':time.process_time()-start,'python':sys.version,'platform':platform.platform(),
    'script_sha256':sha(Path(__file__)),'actual_command':[sys.executable,str(Path(__file__).resolve())]+[z for k,v in vars(args).items() for z in ('--'+k,str(v))]}
args.output.write_text(json.dumps(json_conv(result),indent=2)+'\n')
print(json.dumps({'status':result['status'],'runs':len(summary),'total_main_cpu_seconds':result['total_main_cpu_seconds'],'output_sha256':sha(args.output)}),flush=True)
