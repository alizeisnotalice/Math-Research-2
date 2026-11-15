"""R47 candidate/upstream-domain/outer-y/refcount audit, without min-tree replay.
Five new bracket-refine slab catalogues are audited using frozen R45 checker.
R39 catalogue is reused by its already accepted hash. No large bbox loop runs.
"""
from pathlib import Path
from fractions import Fraction as F
from math import isqrt
from collections import Counter,defaultdict
from datetime import datetime,timezone
import argparse,hashlib,importlib.util,json,sys,time
import review_round045_enumerators as prior
MODULE='929cf92baca532b92159ff81542390acbff3688386dd5104decd0ee55a8ed524'
MATH='3e391beace80ec39a668eaeb81de218d8ab2d982da39d919f8f3ead749d5a3b7'
SIX='185064c894cc24fb18db72da588b1ab040d50885dc54e3ec1f096663f780d4a3'
FIVE='f6863ff64889f91b261cff366af33b72927e95505541620c7bfc9be123a406b2'
sha,read,digest=prior.sha,prior.read,prior.digest

def canon(e):
    a,v,sgn,q=F(e[0]),F(e[1]),int(e[2]),F(e[3]);assert q>=0 and sgn in(-1,0,1)
    if not q or not sgn:return(a,v,0,F(0))
    n,d=isqrt(q.numerator),isqrt(q.denominator)
    return(a+sgn*F(n,d),v,0,F(0))if n*n==q.numerator and d*d==q.denominator else(a,v,sgn,q)

def verify_outer(boxes,domain,H,empty,out):
    fl,fr,fb,ft=map(F,domain);bits=out['bits'];scale=1<<bits;cache={}
    def interval(e):
        a,v,sgn,q=e
        if not sgn:return a,a
        if q not in cache:
            k=isqrt(q.numerator*scale*scale//q.denominator);lo=F(k,scale);hi=lo if lo*lo==q else F(k+1,scale);assert lo*lo<=q<=hi*hi;cache[q]=lo,hi
        lo,hi=cache[q];return(a+lo,a+hi)if sgn==1 else(a-hi,a-lo)
    left,right=canon([fl,0,0,0]),canon([fr,0,0,0]);allx={left,right};family=[];supports=[];permanent=[]
    for j,box in enumerate(boxes):
        l,r,y0,y1=map(canon,box);assert l[1]in(0,1)and r[1]in(0,1)and y0[1]==y1[1]==0;family.append((l,r,y0,y1));allx.update((l,r));lo=max(fb,interval(y0)[0]);hi=min(ft,interval(y1)[1])
        if lo<hi:supports.append((j,lo,hi))
        else:permanent.append(j)
    actual=[(int(j),F(lo),F(hi))for j,lo,hi in out['outer_y_support_catalog']];assert supports==actual and [z['box_index']for z in out['permanently_y_empty_boxes']]==permanent
    for z in out['permanently_y_empty_boxes']:
        j=z['box_index'];assert list(map(F,z['outer_clipped_support']))==[max(fb,interval(family[j][2])[0]),min(ft,interval(family[j][3])[1])]
    order=list(map(canon,out['endpoint_order']));assert len(order)==len(allx)and len(set(order))==len(order)and set(order)==allx
    coverage=out['coverage'];assert coverage['candidate_boxes']==len(family)and coverage['positive_outer_y_supports']==len(supports)and coverage['permanently_y_empty']==len(permanent)and coverage['all_raw_x_endpoint_classes']==len(order)and out['input_candidate_boxes_sha256']==digest(boxes)
    events=defaultdict(list);events[fb];events[ft]
    for j,lo,hi in supports:events[lo].append((j,1));events[hi].append((j,-1))
    levels=sorted(events);counts=Counter({left:1,right:1});active=set();max_ref=1;target=out['limiting_pair'];target_seen=False
    for k,y in enumerate(levels):
        for j,delta in events[y]:
            if delta>0:assert j not in active;active.add(j)
            else:assert j in active;active.remove(j)
            for e in family[j][:2]:counts[e]+=delta;assert counts[e]>=0;max_ref=max(max_ref,counts[e])
        assert counts[left]>=1 and counts[right]>=1
        if k+1==len(levels):continue
        upper=levels[k+1];assert y<upper
        if target is not None and[y,upper]==list(map(F,target['shared_outer_y_strip'])):
            target_seen=True;assert len(active)==target['active_box_count_in_strip'];i,j=target['endpoint_ranks'];assert[order[i],order[j]]==list(map(canon,target['endpoints']))and abs(order[i][1]-order[j][1])==target['relative_speed']==1
            for e,(owner,side)in zip((order[i],order[j]),target['endpoint_owners']):
                if owner=='domain':assert e==({'left':left,'right':right}[side])
                else:assert owner in active and e==family[owner][0 if side=='left'else 1]
            gap=interval(order[j])[0]-interval(order[i])[1];assert gap==F(target['gap_lower'])==F(out['overlap_gap_lower'])>0
    assert not active and {e:n for e,n in counts.items()if n}=={left:1,right:1}and(target_seen or target is None)
    assert coverage['y_levels']==len(levels)and coverage['positive_height_y_strips']==len(levels)-1 and coverage['box_enter_exit_events']==2*len(supports)and coverage['canonical_endpoint_count_updates']==4*len(supports)+2
    assert all(coverage[k]is True for k in('same_y_events_grouped','all_candidate_raw_endpoints_retained','Fx_boundaries_active_throughout_Fy','canonical_multiplicity_restored'))and coverage['velocity_domain_validated']==[0,1]
    H,empty=F(H),F(empty);assert F(out['H'])==H and F(out['delta_empty'])==empty;radius=None if out['overlap_radius']is None else F(out['overlap_radius']);constraints=[('remaining_source_H',H),('remaining_original_empty_radius',empty)]
    if radius is not None:assert radius==F(out['overlap_gap_lower'])/2;constraints.append(('overlap_coupled_endpoint_order',radius))
    assert out['radius_constraints']==[[k,str(v)]for k,v in constraints];reason,final=min(constraints,key=lambda z:z[1]);assert F(out['final_delta'])==final and out['final_delta_limiter']==reason
    oldraw=None if out['old_global_raw_radius']is None else F(out['old_global_raw_radius']);oldfinal=min([H,empty]+([]if oldraw is None else[oldraw]));assert F(out['old_global_final_delta'])==oldfinal and F(out['improvement_factor'])==final/oldfinal>=1
    return {'candidate_boxes':len(family),'outer_supports_independently_reconstructed':len(supports),'permanently_y_empty':len(permanent),'canonical_x_classes':len(order),'positive_y_strips':len(levels)-1,'reference_updates_verified':4*len(supports)+2,'max_shared_reference_count':max_ref,'Fx_remains_active_and_references_balance':True,'limiter_owner_strip_verified':target_seen,'final_delta':str(final),'improvement_factor':str(final/oldfinal),'limiter':reason,'global_minimum_replayed_by_this_script':False}

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--six-index',type=Path,required=True);ap.add_argument('--five-index',type=Path,required=True);ap.add_argument('--topic',type=Path,required=True);ap.add_argument('--geometry-topic',type=Path,required=True);ap.add_argument('--module',type=Path,required=True);ap.add_argument('--math',type=Path,required=True);ap.add_argument('--sweep',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();start=time.perf_counter();assert sha(a.six_index)==SIX and sha(a.five_index)==FIVE and sha(a.module)==MODULE and sha(a.math)==MATH and sha(Path(prior.__file__))=='5075bf81a2c0c0c137e1c305395c2254f047943c415ee0d495f84bd6353a2279';six=read(a.six_index)['records'];five=read(a.five_index)['runs'];assert len(six)==6 and len(five)==5
    basep=a.topic/'computation/inputs/r041_final_0036.json';base=read(basep);assert sha(basep)==prior.BASE;b8p=a.topic/'computation/inputs/r008_final_0611.json';P8=F(read(b8p)['Px']);new_catalogues=[];centers={'q0':F(-139,172800000),'q1':F(-21,57600000),'q2':F(-139,172800000),'q3':F(-139,172800000),'all':F(-45,57600000)}
    for r in five:
        name=r['name'];path=Path(r['source_dir']);d=read(path/'slab.json');m=read(path/'run_manifest.json');vel=read(path/'velocities.json');assert all(sha(path/k)==h for k,h in r['hashes'].items());q=None if name=='all'else int(name[-1]);labels=[51*k+j for k in(range(4)if q is None else[q])for j in(26,27)];assert vel==[int(j in labels)for j in range(204)]and d['direction_file_sha256']==sha(path/'velocities.json')==m['velocity_sha256']
        tau=centers[name];assert F(r['normalized_s0'])==F(m['normalized_s0'])==tau and F(m['physical_s0'])==F(d['s0'])==tau*P8 and F(m['normalized_H'])==F(1,100000000)and F(m['physical_H'])==F(d['requested_H'])==P8/100000000 and m['input_sha256']==prior.BASE
        runner=a.geometry_topic/'computation'/('run_r45_bracket_'+name+'.py');assert sha(runner)==m['runner_sha256'];check=prior.check_slab(d,base,vel);check.update({'name':name,'slab_sha256':sha(path/'slab.json'),'normalized_s0':str(tau),'physical_s0':d['s0'],'labels':labels,'runner_sha256':sha(runner)});assert r['counts']==d['counts'];new_catalogues.append(check);print(name,'new catalogue passed',flush=True)
    records=[]
    for r in six:
        p=Path(r['source_slab']);o=Path(r['output']);assert sha(p)==r['source_slab_sha256']and sha(o)==r['output_sha256'];d=read(p);out=read(o);assert out['code_sha256']==r['code_sha256']==MODULE and out['mathematical_contract_sha256']==MATH and out['input_file_sha256']==sha(p)and out['sweep_module_sha256']==sha(a.sweep)and out['interval_primitives_sha256']==sha(a.sweep.with_name('r039_interval_primitives.py'))
        if r['name']=='r039':assert sha(p)=='c0d262c39a01db8d50cd2fee1a00c7867aa8e08417c4fc85522103b70ce13b25'
        else:assert any(z['name']==r['name']and z['slab_sha256']==sha(p)for z in new_catalogues)
        assert d['schema']=='position-slab-enumerator-v1'and out['source_metadata']=={k:d[k]for k in('schema','s0','code_sha256','source_catalog_sha256','input_source_catalog_sha256')if k in d};check=verify_outer(d['boxes'],d['domain'],d['H'],d['delta_empty'],out)
        for k in('old_global_final_delta','overlap_radius','final_delta','improvement_factor'):assert out[k]==r[k]
        assert out['final_delta_limiter']==r['final_limiter'];check.update({'name':r['name'],'source_slab_sha256':sha(p),'output_sha256':sha(o)});records.append(check);print(r['name'],'outer support/reference/domain passed',flush=True)
    # Bounded behavioral checks specific to omitted/currently-empty candidates.
    spec=importlib.util.spec_from_file_location('r47_reviewed_overlap',a.module);module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    def e(x,v=0):return[str(x),str(v),0,'0']
    def box(l,r,y,z,vl=0,vr=0):return[e(l,vl),e(r,vr),e(y),e(z)]
    tests=[]
    families=[('current_x_empty_kept',[box(2,1,0,1,0,1)],[-10,10,-4,4],F(3),F(3),F(1,2)),('shared_node_after_exit',[box(0,2,0,1),box(0,2,0,2),box(F(201,100),F(301,100),1,2,1,1)],[-10,10,-4,4],F(3),F(3),F(1,200)),('Fx_shared_node_stays',[box(0,1,0,1),box(F(999,100),F(1099,100),2,3,1,1)],[0,10,-4,4],F(3),F(3),F(1,200)),('remaining_original_domain_caps',[box(2,1,0,1,0,1)],[-10,10,-4,4],F(1,40),F(1,200),F(1,200)),('fixed_y_empty_accounted',[box(0,1,9,10)],[-10,10,-4,4],F(3),F(2),F(2))]
    for name,boxes,domain,H,empty,expected in families:
        out=module.overlap_radius(boxes,domain,H,empty,sweep_module=a.sweep,bits=80);check=verify_outer(boxes,domain,H,empty,out);assert F(out['final_delta'])==expected;check['name']=name;tests.append(check)
    result={'status':'passed','scope':'candidate and upstream-domain/outer-y/event/refcount/Fx audit; six minimum-gap trees and beta certificates are independently reviewed by envelope','math_sha256':MATH,'module_sha256':MODULE,'six_index_sha256':SIX,'five_new_catalogue_index_sha256':FIVE,'R45_checker_reused_sha256':sha(Path(prior.__file__)),'R39_catalogue_reused_sha256':'c0d262c39a01db8d50cd2fee1a00c7867aa8e08417c4fc85522103b70ce13b25','five_new_catalogues':new_catalogues,'six_support_attachments':records,'bounded_behavior_checks':tests,'total_new_active_boxes_verified':sum(z['active_boxes_verified']for z in new_catalogues),'large_bbox_loops_replayed':False,'recenter_remaining_domain_automatically_validated':False,'recenter_condition':'caller must independently establish remaining_source_H and remaining_original_empty_radius for the full original candidate family; module only preserves the supplied caps','script_sha256':sha(Path(__file__)),'python':sys.version,'command':[sys.executable,*sys.argv],'verified_utc':datetime.now(timezone.utc).isoformat(),'elapsed_seconds':time.perf_counter()-start}
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items()if k not in('five_new_catalogues','six_support_attachments','bounded_behavior_checks')},indent=2))
if __name__=='__main__':main()
