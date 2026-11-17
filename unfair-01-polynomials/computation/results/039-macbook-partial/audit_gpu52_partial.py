#!/usr/bin/env python3
"""Independent audit of an interrupted GPU52 stream, never a completion claim."""
from __future__ import annotations
import argparse, base64, copy, hashlib, json, os, sys, time
from pathlib import Path
import numpy as np

CONTROL_A={0,1,3}; CONTROL_B={0,5,11,20,49}
def sha_bytes(b): return hashlib.sha256(b).hexdigest()
def sha(p): return sha_bytes(p.read_bytes())
def support(a,tol=1e-12): return {i for i,x in enumerate(a) if abs(x)>tol}

def group_poly(group):
    roots=np.asarray([complex(*z) for z in group],dtype=np.complex128)
    p=np.polynomial.polynomial.polyfromroots(roots); p=p/p[0]
    assert np.max(abs(p.imag))<1e-6
    return p.real

def factor_pair(task,mask):
    groups=task['groups']; assert 0<=mask<task['masks']==2**(groups-1)
    bits=[0]+list(map(int,format(mask,f'0{groups-1}b')))
    a=np.array([1.]); b=np.array([1.])
    for bit,group in zip(bits,task['roots']):
        q=group_poly(group)
        if bit:a=np.polynomial.polynomial.polymul(a,q)
        else:b=np.polynomial.polynomial.polymul(b,q)
    return np.pad(a,(0,53-len(a))),np.pad(b,(0,53-len(b)))

def expected_inputs(data):
    assert data['seed']==2026111616 and data['L']==64 and data['batch']==32768 and len(data['tasks'])==24
    rng=np.random.default_rng(2026111616);out=[]
    for task_id in range(24):
        c=rng.integers(0,2,53); c[0]=c[-1]=1
        if task_id==0:
            c[:]=0
            for i in CONTROL_A:
                for j in CONTROL_B:c[i+j]=1
        out.append(c.tolist())
    return out

def validate_inputs(data,search):
    expected=expected_inputs(data); theta=2*np.pi*np.outer(np.arange(64),np.arange(53))/64
    assert np.max(abs(np.asarray(data['DC'])-np.cos(theta)/64/(1.031**np.arange(53))))<1e-14
    assert np.max(abs(np.asarray(data['DS'])-np.sin(theta)/64/(1.031**np.arange(53))))<1e-14
    assert all(x in search for x in ['GPUArrayQ','gpu_fallback_rejected','Range[offset,Min[total-1','IntegerDigits[#,2,ng-1]'])
    root_error=log_error=0.
    for tid,task in enumerate(data['tasks']):
        assert task['id']==tid and task['C']==expected[tid] and len(task['C'])==53 and set(task['C'])<={0,1}
        assert task['C'][0]==task['C'][-1]==1 and task['groups']==len(task['roots']) and task['masks']==2**(task['groups']-1)
        assert len(task['LR'])==len(task['LI'])==task['groups']
        roots=[complex(*z) for g in task['roots'] for z in g]; assert len(roots)==52
        product=np.polynomial.polynomial.polyfromroots(roots); product=product/product[0]
        root_error=max(root_error,float(np.max(abs(product-np.asarray(task['C'])))))
        points=1.031*np.exp(2j*np.pi*np.arange(64)/64)
        for group,lr,li in zip(task['roots'],task['LR'],task['LI']):
            v=np.ones(64,dtype=complex)
            for z in group:v*=(points-complex(*z))/(-complex(*z))
            assert np.min(abs(v))>1e-10
            got=np.log(v); log_error=max(log_error,float(max(np.max(abs(got.real-np.asarray(lr))),np.max(abs(got.imag-np.asarray(li))))))
    assert root_error<3e-5 and log_error<1e-11
    return {'max_root_product_error':root_error,'max_log_input_error':log_error}

def parse_events(path):
    rows=[]
    for raw in path.read_text().splitlines():
        if raw.startswith('{'): rows.append(json.loads(raw))
    assert rows and rows[0]['event']=='controller_started'
    return rows

def expected_kernel_hash(data,search,task):
    one={**data,'tasks':[task]}
    wl='data=ImportString['+json.dumps(json.dumps(one,separators=(',',':')))+',"RawJSON"];\n'+search
    return sha_bytes(wl.encode()),len(wl.encode())

def validate_payload(inputs,search,remote,payload):
    rebuilt='INPUT_B64='+repr(base64.b64encode(inputs).decode())+'\nSEARCH_B64='+repr(base64.b64encode(search.encode()).decode())+'\n'+remote
    assert payload==rebuilt
    return sha_bytes(rebuilt.encode())

def audit_events(data,rows,search):
    assert rows[0]['tasks']==24 and rows[0]['seconds_limit']==14400 and rows[0]['storage']=='RAM only; no remote task files'
    byid={t['id']:t for t in data['tasks']}; cursors={}; samples=count=0; max_coef=max_prod=0.; near=[]; control=[]
    current=None; started=[]; completed=[]; reaped=[]
    for e in rows[1:]:
        kind=e['event']
        if kind=='kernel_started':
            assert current is None and e['task']==len(started) and e['task'] in byid
            wanted_h,wanted_n=expected_kernel_hash(data,search,byid[e['task']])
            assert e['program_sha256']==wanted_h and e['program_bytes']==wanted_n
            current=e['task']; cursors[current]=0; started.append(current)
        elif kind=='started':
            assert current is not None and e['tasks']==1 and e['batch']==32768 and e['pid']
        elif kind=='batch':
            assert current is not None and e['task']==current
            task=byid[current]; assert e['first']==cursors[current]
            assert e['count']==e['last']-e['first']+1 and 1<=e['count']<=32768 and e['last']<task['masks'] and e['gpu_arrays_verified'] is True
            cursors[current]=e['last']+1; count+=e['count']; assert 1<=len(e['cases'])<=4
            for case in e['cases']:
                assert set(case)=={'mask','score','A','B'} and e['first']<=case['mask']<=e['last'] and case['mask'] and np.isfinite(case['score']) and case['score']>=0
                a,b=factor_pair(task,case['mask']); sa=np.asarray(case['A']); sb=np.asarray(case['B']); assert sa.shape==sb.shape==(53,)
                max_coef=max(max_coef,float(max(np.max(abs(a-sa)),np.max(abs(b-sb)))))
                prod=np.polynomial.polynomial.polymul(a,b); max_prod=max(max_prod,float(np.max(abs(prod-np.asarray(task['C'])))))
                score=max(0.,-float(np.min(a)),-float(np.min(b)),float(np.max(a))-1.,float(np.max(b))-1.)
                if score<2e-5:near.append({'task':current,'mask':case['mask'],'cpu_score':score,'gpu_score':float(case['score'])})
                rounded_a=np.rint(sa); rounded_b=np.rint(sb)
                rounding=max(np.max(abs(sa-rounded_a)),np.max(abs(sb-rounded_b)))
                if current==0 and rounding<2e-5 and support(rounded_a)==CONTROL_A and support(rounded_b)==CONTROL_B: control.append((case,sa,sb))
                samples+=1
        elif kind=='complete':
            assert current is not None and cursors[current]==byid[current]['masks'] and e['seconds']>0
            completed.append(current)
        elif kind=='kernel_reaped':
            assert current is not None and current in completed and e['task']==current and e['returncode']==0 and e['timed_out'] is False
            reaped.append(current); current=None
        else: raise AssertionError('unexpected event '+kind)
    assert started==[0,1,2] and completed==[0,1] and reaped==[0,1] and current==2 and cursors[2]<byid[2]['masks']
    assert max_coef<3e-4 and max_prod<2e-4 and len(control)==1, (max_coef,max_prod,len(control))
    case,sa,sb=control[0]; ra=np.rint(sa); rb=np.rint(sb)
    assert support(ra)==CONTROL_A and support(rb)==CONTROL_B and np.max(abs(np.polynomial.polynomial.polymul(ra,rb)-np.asarray(byid[0]['C'])))==0
    return {'received_tasks_complete':completed,'received_task_prefix':{'task':2,'masks':cursors[2],'of':byid[2]['masks']},'received_masks':count,'saved_cases':samples,'max_gpu_cpu_coefficient_error':max_coef,'max_sample_product_error':max_prod,'near_saved_cpu_score_lt_2e-5':near,'fair_control':{'mask':case['mask'],'exact_rounded_product_binary':True,'rounding_change':float(max(np.max(abs(sa-ra)),np.max(abs(sb-rb))))}}

def rejected(label,thunk):
    try: thunk()
    except (AssertionError,KeyError,IndexError,TypeError): return label
    raise AssertionError('audit accepted mutation '+label)

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--root',type=Path,required=True); ap.add_argument('--out',type=Path,required=True); a=ap.parse_args(); assert os.path.ismount('/Volumes/data'); started=time.monotonic()
    root=a.root; inputs=(root/'inputs.json').read_bytes(); data=json.loads(inputs); search=(root/'search.wl').read_text(); remote=(root/'remote_gpu.py').read_text(); payload=(root/'payload.py').read_text(); rows=parse_events(root/'stdout.jsonl')
    inp=validate_inputs(data,search); payload_hash=validate_payload(inputs,search,remote,payload); stream=audit_events(data,rows,search)
    bad_input=copy.deepcopy(data); bad_input['tasks'][0]['C'][0]=0
    bad_masks=copy.deepcopy(data); bad_masks['tasks'][0]['masks']-=1
    bad_jump=copy.deepcopy(rows); next(x for x in bad_jump if x['event']=='batch')['first']=1
    bad_drop=copy.deepcopy(rows); bad_drop.pop(next(i for i,x in enumerate(bad_drop) if x['event']=='batch' and x['task']==2))
    bad_case=copy.deepcopy(rows); next(x for x in bad_case if x['event']=='batch')['cases'][0]['A'][0]+=.1
    bad_complete=copy.deepcopy(rows); bad_complete.append({'event':'complete','seconds':1.})
    negatives=[rejected('input_endpoint',lambda:validate_inputs(bad_input,search)),rejected('wrong_mask_total',lambda:audit_events(bad_masks,rows,search)),rejected('batch_cursor_jump',lambda:audit_events(data,bad_jump,search)),rejected('dropped_partial_batch',lambda:audit_events(data,bad_drop,search)),rejected('saved_coefficient',lambda:audit_events(data,bad_case,search)),rejected('fake_completion',lambda:audit_events(data,bad_complete,search))]
    out={'status':'passed_interrupted_gpu52_partial_audit','scope':'Only two whole task mask spaces and the recorded contiguous task-2 prefix are audited. The global batch was interrupted; all unreceived masks are unknown. Numerical values are not exact exclusions or factorization certificates.','input':inp,'stream':stream,'payload_sha256':payload_hash,'negative_tests_rejected':negatives,'source_sha256':{'inputs':sha(root/'inputs.json'),'search':sha(root/'search.wl'),'remote':sha(root/'remote_gpu.py'),'payload':sha(root/'payload.py'),'stdout':sha(root/'stdout.jsonl'),'auditor':sha(Path(__file__))},'interpreter':sys.version,'numpy':np.__version__,'elapsed_seconds':time.monotonic()-started}
    a.out.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n'); print(json.dumps({'status':out['status'],'stream':stream,'negative_tests_rejected':negatives,'seconds':out['elapsed_seconds']},sort_keys=True))
if __name__=='__main__':main()
