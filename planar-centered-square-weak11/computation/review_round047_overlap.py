#!/usr/bin/env python3
"""R47 independent overlap-gap heap audit plus 5 genuinely new slope replays.
No production CoupledGapTree or candidate-bbox enumeration is used for real data.
"""
import argparse,ast,bisect,gzip,hashlib,heapq,importlib.util,json,platform,sys,time
from fractions import Fraction as F
from functools import lru_cache
from math import isqrt,lcm
from pathlib import Path
from types import SimpleNamespace
sys.dont_write_bytecode=True
parser=argparse.ArgumentParser(description=__doc__)
for key in ['index','new-index','module','contract','sweep-module','old-kernel','old-replay','prior-check','output']:
    parser.add_argument('--'+key,required=True,type=Path)
args=parser.parse_args()
for k,v in vars(args).items():setattr(args,k,v.resolve())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def serial(x):
    if isinstance(x,F):return str(x)
    if isinstance(x,dict):return {str(k):serial(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):return [serial(v) for v in x]
    return x
pins={args.index:'185064c894cc24fb18db72da588b1ab040d50885dc54e3ec1f096663f780d4a3',
 args.new_index:'f6863ff64889f91b261cff366af33b72927e95505541620c7bfc9be123a406b2',
 args.module:'929cf92baca532b92159ff81542390acbff3688386dd5104decd0ee55a8ed524',
 args.contract:'3e391beace80ec39a668eaeb81de218d8ab2d982da39d919f8f3ead749d5a3b7',
 args.old_kernel:'46ac35f415025e03af461b9ea26a71f1c2a71c04968389ef1dbc66a6f4057bd2',
 args.old_replay:'6f296cca49b3cd6c6bbdd62e837ce2d678727424cd8fadc8e613db985441ff19',
 args.prior_check:'2667aa8d1803bfe51ff11bfaddd9a735405da51ac2a1c1a508b6be511ecd7219'}
for p,h in pins.items():assert sha(p)==h
prior=json.loads(args.prior_check.read_text());P=F(prior['real']['normalization']['base_Px'])
source_sha={name:sha(args.sweep_module.with_name(name)) for name in prior['source_sha256']}
assert source_sha==prior['source_sha256']
def extract(path,names,namespace):
    nodes=[n for n in ast.parse(path.read_text()).body if isinstance(n,ast.FunctionDef) and n.name in names]
    assert len(nodes)==len(names)
    exec(compile(ast.Module(body=nodes,type_ignores=[]),str(path),'exec'),namespace)
ns={'F':F,'isqrt':isqrt,'lru_cache':lru_cache}
extract(args.old_kernel,['root','canon','enclosure','merge'],ns)
canon,enclosure,merge=(ns[k] for k in ['canon','enclosure','merge'])
ns.update({'time':time,'json':json,'gzip':gzip,'hashlib':hashlib,'sha':sha,'Path':Path,'lcm':lcm,'P':P,'source_sha':source_sha})
extract(args.old_replay,['loadraw','replay'],ns)
start=time.process_time();index=json.loads(args.index.read_text());new_index=json.loads(args.new_index.read_text())
assert len(index['records'])==6 and len(new_index['runs'])==5
new_slopes=[]
for entry in new_index['runs']:
    folder=Path(entry['source_dir']);manifest=json.loads((folder/'run_manifest.json').read_text());slope=json.loads((folder/'slope.json').read_text())
    adapted=dict(entry,name=folder.name,normalized_requested_H=manifest['normalized_H'],tail=slope['tail'])
    ns['args']=SimpleNamespace(index=folder.parent/'bracket_refine_index.json')
    item=ns['replay'](adapted);item['direction']=entry['name'];new_slopes.append(item)

# Independent ordered-list/linked-neighbour/heap implementation; no segment tree.
def heap_replay(boxes,domain,bits,supports):
    fl,fr,fb,ft=map(F,domain);left,right=canon([fl,0,0,0]),canon([fr,0,0,0])
    family=[tuple(canon(e) for e in box) for box in boxes];allx={left,right}
    for l,r,y0,y1 in family:
        assert l[1] in [0,1] and r[1] in [0,1] and not y0[1] and not y1[1];allx.update((l,r))
    rows=sorted((enclosure(e,bits)[0],enclosure(e,bits)[1],e) for e in allx)
    assert all(a[1]<b[0] for a,b in zip(rows,rows[1:]))
    order=[row[2] for row in rows];rank={e:i for i,e in enumerate(order)};N=len(order)
    count=[0]*N;previous=[None]*N;following=[None]*N;ordered=[];heap=[];gaps={}
    def gap(i,j):
        key=i,j
        if key not in gaps:gaps[key]=rows[j][0]-rows[i][1];assert gaps[key]>0
        return gaps[key]
    def seam(i,j):
        if i is not None and j is not None and order[i][1]!=order[j][1]:heapq.heappush(heap,(gap(i,j),i,j))
    def update(i,delta):
        old=count[i];count[i]+=delta;assert count[i]>=0
        if bool(old)==bool(count[i]):return
        pos=bisect.bisect_left(ordered,i)
        if count[i]:
            a=ordered[pos-1] if pos else None;b=ordered[pos] if pos<len(ordered) else None
            ordered.insert(pos,i);previous[i],following[i]=a,b
            if a is not None:following[a]=i
            if b is not None:previous[b]=i
            seam(a,i);seam(i,b)
        else:
            assert ordered[pos]==i;ordered.pop(pos);a,b=previous[i],following[i]
            if a is not None:following[a]=b
            if b is not None:previous[b]=a
            previous[i]=following[i]=None;seam(a,b)
    def minimum():
        while heap:
            g,i,j=heap[0]
            if count[i] and count[j] and following[i]==j:return g,i,j
            heapq.heappop(heap)
        return None
    update(rank[left],1);update(rank[right],1)
    events={fb:[],ft:[]};seen=set()
    for i,a,b in supports:
        assert i not in seen and 0<=i<len(family);seen.add(i);a,b=F(a),F(b);assert fb<=a<b<=ft
        events.setdefault(a,[]).append((i,1));events.setdefault(b,[]).append((i,-1))
    levels=sorted(events);active_boxes=set();best=None;strip=None;changes=0
    for k,y in enumerate(levels):
        for i,d in events[y]:
            if d>0:assert i not in active_boxes;active_boxes.add(i)
            else:assert i in active_boxes;active_boxes.remove(i)
            l,r=family[i][:2];update(rank[l],d);update(rank[r],d);changes+=1
        if k+1<len(levels):
            assert y<levels[k+1];candidate=minimum()
            if candidate is not None and (best is None or candidate<best):best=candidate;strip=(y,levels[k+1])
    assert not active_boxes and changes==2*len(supports)
    expected=[0]*N;expected[rank[left]]+=1;expected[rank[right]]+=1
    assert count==expected
    global_candidates=[((b[0]-a[1])/2,i,i+1) for i,(a,b) in enumerate(zip(rows,rows[1:])) if a[2][1]!=b[2][1]]
    global_best=min(global_candidates) if global_candidates else None
    return {'best':best,'strip':strip,'order':order,'rows':rows,'levels':levels,
            'supports':supports,'family':family,'rank':rank,'global_best':global_best,'events':changes}

# Bounded new tests compare every coactive mixed pair, not the production tree.
spec=importlib.util.spec_from_file_location('audited_overlap_module',args.module);module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
def e(a,v=0,sgn=0,q=0):return [str(a),str(v),sgn,str(q)]
def rect(a,b,c,d,u=0,v=0):return [e(a,u),e(b,v),e(c),e(d)]
domain=['-10','10','-4','4']
toy_cases=[('disjoint',[rect(0,2,0,1),rect('201/100','301/100',2,3,1,1)],5,5),
 ('positive_overlap',[rect(0,2,0,2),rect('201/100','301/100',1,3,1,1)],5,5),
 ('touch_only',[rect(0,2,0,1),rect('201/100','301/100',1,2,1,1)],5,5),
 ('shared_canonical',[rect(0,2,0,1),rect(0,2,0,3),rect('201/100','301/100',2,4,1,1)],5,5),
 ('empty_candidate_preserved',[rect(1,0,0,1,1,0)],1,1),
 ('domain_clip',[rect('999/100','1099/100',0,1,1,1)],1,1),
 ('root_x_disjoint',[[e(0),e(0,0,1,2),e(0),e(1)],[e('1/1000',1,1,2),e(3,1),e(2),e(3)]],5,5),
 ('root_y_false_overlap',[[e(0),e(1),e(0),e(0,0,1,2)],[e('1001/1000',1),e('2001/1000',1),e(0,0,1,2),e(3)]],5,5),
 ('omitted_box_domain_cap',[rect(0,2,0,1),rect('201/100','301/100',2,3,1,1)],F(1,100),F(1,400))]
toys=[]
for name,boxes,H,empty in toy_cases:
    out=module.overlap_radius(boxes,domain,H,empty,sweep_module=args.sweep_module)
    own=heap_replay(boxes,domain,out['bits'],out['outer_y_support_catalog']);best=None
    for y,z in zip(own['levels'],own['levels'][1:]):
        ranks={own['rank'][canon([domain[0],0,0,0])],own['rank'][canon([domain[1],0,0,0])]}
        for i,a,b in own['supports']:
            if F(a)<=y<z<=F(b):ranks.update(own['rank'][x] for x in own['family'][i][:2])
        ranks=sorted(ranks)
        for pos,i in enumerate(ranks):
            for j in ranks[pos+1:]:
                if own['order'][i][1]!=own['order'][j][1]:
                    candidate=(own['rows'][j][0]-own['rows'][i][1],i,j)
                    if best is None or candidate<best:best=candidate
    assert best==own['best']
    radius=best[0]/2 if best else None
    assert radius==(F(out['overlap_radius']) if out['overlap_radius'] is not None else None)
    assert F(out['final_delta'])==min([F(H),F(empty)]+([radius] if radius is not None else []))
    toys.append({'case':name,'overlap_radius':radius,'final_delta':F(out['final_delta']),'all_coactive_pairs_match':True})
for boxes in [[rect(0,1,0,1,0,2)],[rect(0,1,0,1),rect(1,2,2,3,1,1)]]:
    try:module.overlap_radius(boxes,domain,1,1,sweep_module=args.sweep_module)
    except ValueError:pass
    else:raise AssertionError('unsupported input accepted')

real=[]
for record in index['records']:
    t=time.process_time();source_path=Path(record['source_slab']);output_path=Path(record['output'])
    assert sha(source_path)==record['source_slab_sha256'] and sha(output_path)==record['output_sha256']
    data=json.loads(source_path.read_text());out=json.loads(output_path.read_text());assert data['schema']=='position-slab-enumerator-v1'
    assert out['input_file_sha256']==record['source_slab_sha256'] and out['code_sha256']==record['code_sha256']==sha(args.module)
    assert out['mathematical_contract_sha256']==sha(args.contract)
    assert out['sweep_module_sha256']==source_sha['r039_endpoint_sweep.py'] and out['interval_primitives_sha256']==source_sha['r039_interval_primitives.py']
    assert out['input_candidate_boxes_sha256']==hashlib.sha256(json.dumps(data['boxes'],sort_keys=True,separators=(',',':')).encode()).hexdigest()
    own=heap_replay(data['boxes'],data['domain'],out['bits'],out['outer_y_support_catalog'])
    assert own['order']==[canon(e) for e in out['endpoint_order']]
    g,i,j=own['best'];assert g==F(out['overlap_gap_lower']) and [i,j]==out['limiting_pair']['endpoint_ranks']
    assert list(own['strip'])==list(map(F,out['limiting_pair']['shared_outer_y_strip']))
    old=own['global_best'][0];new=g/2;H,empty=F(data['H']),F(data['delta_empty'])
    assert F(out['H'])==H and F(out['delta_empty'])==empty
    old_final=min(old,H,empty);final=min(new,H,empty)
    assert old==F(out['old_global_raw_radius']) and old_final==F(out['old_global_final_delta'])==F(record['old_global_final_delta'])
    assert new==F(out['overlap_radius'])==F(record['overlap_radius'])
    assert final==F(out['final_delta'])==F(record['final_delta']) and final/old_final==F(out['improvement_factor'])==F(record['improvement_factor'])
    assert final>=old_final>0 and new>0
    assert len(own['levels'])-1==out['coverage']['positive_height_y_strips']
    assert own['events']==out['coverage']['box_enter_exit_events']
    slope=(prior['real'] if record['name']=='r039' else next(r for r in new_slopes if r['direction']==record['name']))
    assert F(slope.get('delta',slope.get('delta_physical')))==old_final
    item={'name':record['name'],'source_slab_sha256':record['source_slab_sha256'],'output_sha256':record['output_sha256'],
       'boxes':len(data['boxes']),'x_classes':len(own['order']),'positive_outer_y_strips':len(own['levels'])-1,
       'overlap_gap_lower':g,'overlap_radius':new,'old_final_delta':old_final,'final_delta':final,
       'H':H,'delta_empty':empty,'improvement_factor':final/old_final,'normalized_final_delta':final/P,
       'beta_lower':F(slope['beta_lower']),'beta_upper':F(slope['beta_upper']),
       'new_beta_replay':record['name']!='r039','limiting_pair_ranks':[i,j],'shared_outer_strip':own['strip'],
       'independent_method':'sorted active rank list, linked neighbours, lazy minimum heap; no production tree',
       'cpu_seconds':time.process_time()-t}
    real.append(item);print(json.dumps({'graph':record['name'],'improved':final>old_final,'strips':item['positive_outer_y_strips'],'cpu_seconds':item['cpu_seconds']}),flush=True)
assert {r['name'] for r in real if r['improvement_factor']>1}=={'r039','q0'}
result={'status':'passed_gap_and_slope_replay_pending_catalogue_binding','scope':'new overlap graph gaps and five new center slopes; upstream candidate/outer-Y/omitted-domain proof supplied by separate audit',
    'index_sha256':sha(args.index),'new_center_index_sha256':sha(args.new_index),'contract_sha256':sha(args.contract),'module_sha256':sha(args.module),
    'source_sha256':source_sha,'old_independent_kernel_sha256':sha(args.old_kernel),'old_independent_slope_replay_sha256':sha(args.old_replay),
    'base_horizontal_period':P,'toy_all_pairs':toys,'nonbinary_and_nonpersistent_collision_rejected':True,
    'new_center_slope_replays':new_slopes,'real_radius_replays':real,
    'candidate_bbox_enumerations_rerun':False,'R39_R45_frozen_outputs_modified':False,
    'main_cpu_seconds':time.process_time()-start,'python':sys.version,'platform':platform.platform(),'script_sha256':sha(Path(__file__)),
    'actual_command':[sys.executable,str(Path(__file__).resolve())]+[z for k,v in vars(args).items() for z in ('--'+k,str(v))]}
args.output.write_text(json.dumps(serial(result),indent=2)+'\n');print(json.dumps({'status':result['status'],'main_cpu_seconds':result['main_cpu_seconds'],'output_sha256':sha(args.output)}),flush=True)
