#!/usr/bin/env python3
"""Independent numerical audit for the completed MacBook GPU48 batch.

All input and saved outputs are read locally.  This validates coverage and
float64 reconstruction, but deliberately does not treat floating-point search
as an exact factorization certificate or an exclusion theorem.
"""
from __future__ import annotations
import argparse, copy, hashlib, json, os, sys, time
from pathlib import Path
import numpy as np

CONTROL_MASK=131584
CONTROL_A={0,1,3}; CONTROL_B={0,5,11,20,45}
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def support(v,tol=1e-12):return {i for i,x in enumerate(v) if abs(x)>tol}
def poly_for_group(group):
    roots=np.asarray([complex(*z) for z in group],dtype=np.complex128)
    p=np.polynomial.polynomial.polyfromroots(roots);p=p/p[0]
    assert np.max(abs(p.imag))<1e-7
    return p.real
def factor_pair(task,mask):
    groups=task['groups'];assert 0<=mask<task['masks']==2**(groups-1)
    bits=[0]+list(map(int,format(mask,f'0{groups-1}b')))
    aa=np.array([1.]);bb=np.array([1.])
    for bit,group in zip(bits,task['roots']):
        p=poly_for_group(group)
        if bit:aa=np.polynomial.polynomial.polymul(aa,p)
        else:bb=np.polynomial.polynomial.polymul(bb,p)
    return np.pad(aa,(0,49-len(aa))),np.pad(bb,(0,49-len(bb)))
def expected_inputs(data):
    assert data['seed']==2026111614 and data['L']==64 and data['batch']==32768 and len(data['tasks'])==24
    rng=np.random.default_rng(2026111614);expected=[]
    for ti in range(24):
        c=rng.integers(0,2,49);c[0]=c[-1]=1
        if ti==0:
            c[:]=0
            for i in CONTROL_A:
                for j in CONTROL_B:c[i+j]=1
        expected.append(c.tolist())
    return expected
def validate_inputs(data, search):
    expected=expected_inputs(data);theta=2*np.pi*np.outer(np.arange(64),np.arange(49))/64
    assert np.max(abs(np.asarray(data['DC'])-np.cos(theta)/64/(1.031**np.arange(49))))<1e-14
    assert np.max(abs(np.asarray(data['DS'])-np.sin(theta)/64/(1.031**np.arange(49))))<1e-14
    assert 'GPUArrayQ' in search and 'gpu_fallback_rejected' in search and 'Range[offset,Min[total-1' in search and 'IntegerDigits[#,2,ng-1]' in search
    max_root=0.;max_log=0.
    for ti,task in enumerate(data['tasks']):
        assert task['id']==ti and task['C']==expected[ti] and len(task['C'])==49 and task['C'][0]==task['C'][-1]==1 and set(task['C'])<={0,1}
        assert task['groups']==len(task['roots']) and task['masks']==2**(task['groups']-1) and len(task['LR'])==len(task['LI'])==task['groups']
        allroots=[complex(*z) for group in task['roots'] for z in group];assert len(allroots)==48
        product=np.polynomial.polynomial.polyfromroots(allroots);product=product/product[0]
        err=float(np.max(abs(product-np.asarray(task['C']))));max_root=max(max_root,err);assert err<1e-5
        points=1.031*np.exp(2j*np.pi*np.arange(64)/64)
        for group,lr,li in zip(task['roots'],task['LR'],task['LI']):
            vals=np.ones(64,dtype=np.complex128)
            for z in group:vals*=(points-complex(*z))/(-complex(*z))
            assert np.min(abs(vals))>1e-10
            got=np.log(vals);err=float(max(np.max(abs(got.real-np.asarray(lr))),np.max(abs(got.imag-np.asarray(li)))));max_log=max(max_log,err);assert err<1e-12
    return {'max_root_product_error':max_root,'max_log_input_error':max_log}
def read_events(path):
    events=[]
    for line in path.read_text().splitlines():
        if line.startswith('{'):events.append(json.loads(line))
    assert events and events[0].get('event')=='started' and events[-1].get('event')=='complete'
    return events
def validate_execution_binding(input_text,search,exit_record):
    """Rebuild exactly the Wolfram stdin payload used by the documented bridge."""
    payload='data=ImportString['+json.dumps(input_text)+',"RawJSON"];\n'+search
    assert exit_record['returncode']==0
    assert exit_record['input_sha256']==hashlib.sha256(input_text.encode()).hexdigest()
    assert exit_record['program_sha256']==hashlib.sha256(payload.encode()).hexdigest()
    assert exit_record['command'][0]=='ssh'
    return hashlib.sha256(payload.encode()).hexdigest()
def validate_events(data,events):
    assert events[0]=={'event':'started','pid':events[0]['pid'],'version':'15.0.0 for Mac OS X ARM (64-bit) (May 15, 2026)','tasks':24,'batch':32768}
    assert events[-1]['event']=='complete' and events[-1]['pid']==events[0]['pid'] and events[-1]['seconds']>0
    cursors={t['id']:0 for t in data['tasks']};count=samples=0;max_gpu_cpu=0.;max_product=0.;best=float('inf');best_noncontrol=float('inf');near=[];control=[]
    for e in events[1:-1]:
        assert set(e)=={'event','task','first','last','count','gpu_arrays_verified','cases','seconds'} and e['event']=='batch'
        task=data['tasks'][e['task']];assert e['first']==cursors[e['task']] and e['count']==e['last']-e['first']+1 and 1<=e['count']<=data['batch'] and e['last']<task['masks'] and e['gpu_arrays_verified'] is True and e['seconds']>=0
        cursors[e['task']]=e['last']+1;count+=e['count'];assert isinstance(e['cases'],list) and 1<=len(e['cases'])<=4
        for case in e['cases']:
            assert set(case)=={'mask','score','A','B'} and e['first']<=case['mask']<=e['last'] and case['mask']!=0
            assert np.isfinite(case['score']) and case['score']>=0
            aa,bb=factor_pair(task,case['mask']);saved_a=np.asarray(case['A']);saved_b=np.asarray(case['B']);assert saved_a.shape==saved_b.shape==(49,)
            err=float(max(np.max(abs(aa-saved_a)),np.max(abs(bb-saved_b))));max_gpu_cpu=max(max_gpu_cpu,err);assert err<2e-4
            # Degree-48 roots reconstructed from the supplied numerical root
            # lists incur conditioning error; this is an audit tolerance, not
            # an exact polynomial identity claim.  Its value is still far
            # below the GPU sample tolerance above.
            product=np.polynomial.polynomial.polymul(aa,bb);cerr=float(np.max(abs(product-np.asarray(task['C']))));max_product=max(max_product,cerr);assert cerr<2e-5
            score=max(0.,-float(np.min(aa)),-float(np.min(bb)),float(np.max(aa))-1.,float(np.max(bb))-1.)
            best=min(best,score);samples+=1
            if task['id']!=0: best_noncontrol=min(best_noncontrol,score)
            if score<2e-5:
                near.append({'task':task['id'],'mask':case['mask'],'cpu_score':score,'gpu_saved_score':float(case['score']),'fair_control':task['id']==0 and case['mask']==CONTROL_MASK})
            if task['id']==0 and case['mask']==CONTROL_MASK:control.append((aa,bb,saved_a,saved_b,err,score,float(case['score'])))
    assert all(cursors[t['id']]==t['masks'] for t in data['tasks'])
    assert len(control)==1
    aa,bb,sa,sb,err,score,gpu_score=control[0]
    expected_a=np.zeros(49);expected_b=np.zeros(49)
    for i in CONTROL_A:expected_a[i]=1
    for i in CONTROL_B:expected_b[i]=1
    assert support(expected_a)==CONTROL_A and support(expected_b)==CONTROL_B
    assert np.max(abs(np.polynomial.polynomial.polymul(expected_a,expected_b)-np.asarray(data['tasks'][0]['C'])))==0
    control_rounding=float(max(np.max(abs(aa-expected_a)),np.max(abs(bb-expected_b)),np.max(abs(sa-expected_a)),np.max(abs(sb-expected_b))))
    assert control_rounding<2e-5
    assert np.isfinite(best_noncontrol)
    return {'masks':count,'samples':samples,'max_gpu_cpu_coefficient_error':max_gpu_cpu,'max_sample_product_error':max_product,'minimum_cpu_score':best,'minimum_noncontrol_cpu_score':best_noncontrol,'near_saved_cpu_score_lt_2e-5':near,'control':{'mask':CONTROL_MASK,'A_support':sorted(CONTROL_A),'B_support':sorted(CONTROL_B),'exact_product_binary':True,'maximum_rounding_change':control_rounding,'gpu_score':gpu_score,'cpu_score':score}}
def rejected(label,thunk):
    try:thunk()
    except (AssertionError,KeyError,TypeError,IndexError):return label
    raise AssertionError('production audit accepted mutation: '+label)
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--input',type=Path,required=True);ap.add_argument('--stdout',type=Path,required=True);ap.add_argument('--search',type=Path,required=True);ap.add_argument('--exit',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);x=ap.parse_args();assert os.path.ismount('/Volumes/data');start=time.monotonic()
    input_text=x.input.read_text();data=json.loads(input_text);events=read_events(x.stdout);search=x.search.read_text();exit_record=json.loads(x.exit.read_text())
    payload_hash=validate_execution_binding(input_text,search,exit_record)
    input_report=validate_inputs(data,search);event_report=validate_events(data,events)
    bad_input=copy.deepcopy(data);bad_input['tasks'][0]['C'][0]=0
    bad_mask=copy.deepcopy(data);bad_mask['tasks'][0]['masks']-=1
    bad_events=copy.deepcopy(events);bad_events[1]['first']=1
    bad_sample=copy.deepcopy(events);bad_sample[1]['cases'][0]['A'][0]+=0.1
    bad_exit=copy.deepcopy(exit_record);bad_exit['program_sha256']='0'*64
    negatives=[rejected('input_nonbinary_endpoint',lambda:validate_inputs(bad_input,search)),rejected('input_wrong_mask_total',lambda:validate_events(bad_mask,events)),rejected('event_noncontinuous_cursor',lambda:validate_events(data,bad_events)),rejected('saved_coefficient_mismatch',lambda:validate_events(data,bad_sample)),rejected('stdin_payload_hash_mismatch',lambda:validate_execution_binding(input_text,search,bad_exit))]
    out={'status':'passed_gpu48_independent_numerical_audit','input':input_report,'events':event_report,'remote_stdin_payload_sha256':payload_hash,'negative_tests_rejected':negatives,'coverage':'every declared mask interval is contiguous and all task cursors terminate exactly at their declared mask counts','scope':'Float64 reconstructions of float32 GPU samples; the fair control is exact separately. No numerical result is an exact unfair factorization, proof, or exclusion.','source_sha256':{'inputs':sha(x.input),'stdout':sha(x.stdout),'search':sha(x.search),'exit':sha(x.exit),'auditor':sha(Path(__file__))},'interpreter':sys.version,'numpy':np.__version__,'elapsed_seconds':time.monotonic()-start}
    x.out.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n');print(json.dumps(out,sort_keys=True))
if __name__=='__main__':main()
