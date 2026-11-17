#!/usr/bin/env python3
"""Single-remote-process Metal scheduler with a durable remote JSONL stream.

This is deliberately a separate handoff candidate, not a live replacement for
GPU v3.  One Swift process loads an input and builds the Metal pipelines once,
then emits one fsynced, bound top-four event per requested batch.  The local
controller validates that stream before writing its own fsynced ledger.
"""
import argparse, hashlib, json, math, os, subprocess, sys, time
from pathlib import Path

HERE=Path(__file__).resolve().parent
TEMP=Path('/Volumes/data/CodexTemp/T2/unfair-01-polynomials')
INPUTS=TEMP/'round038-macbook-gpu52-direct/inputs.json'
SOURCE=HERE/'NativeScoreProduction.swift'
REMOTE='/Users/zhengzhihao/CodexRemoteMath/a30-native-production42b'
SSH=['ssh','-o','HostName=fe80::8c00:aeff:feee:2fe%%anpi1','-o','HostKeyAlias=169.254.63.75','-o','BatchMode=yes','-o','StrictHostKeyChecking=yes','-o','ConnectTimeout=10','macbook-math']
SCP=['scp','-o','HostName=fe80::8c00:aeff:feee:2fe%%anpi1','-o','HostKeyAlias=169.254.63.75','-o','BatchMode=yes','-o','StrictHostKeyChecking=yes','-o','ConnectTimeout=10']

def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def fsync_file(p):
    with open(p,'rb') as f: os.fsync(f.fileno())
def write_json(p,x):
    with open(p,'w') as f: json.dump(x,f,sort_keys=True,indent=2);f.write('\n');f.flush();os.fsync(f.fileno())
def append(p,x):
    with open(p,'a') as f:f.write(json.dumps(x,sort_keys=True,separators=(',',':'))+'\n');f.flush();os.fsync(f.fileno())
def remote(command, capture=True): return subprocess.run(SSH+[command],check=True,text=True,capture_output=capture)
def strict_int(x): return type(x) is int

def cpu_poly(group):
    p=[1+0j]
    for rr,ii in group:
        z=complex(rr,ii); q=[0j]*(len(p)+1)
        for i,c in enumerate(p): q[i]+=-z*c;q[i+1]+=c
        p=q
    p=[x/p[0] for x in p]
    if not all(math.isfinite(x.real) and math.isfinite(x.imag) for x in p) or max(abs(x.imag) for x in p)>=1e-7: raise ValueError('invalid real root group')
    return [x.real for x in p]
def conv(a,b):
    o=[0.]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):o[i+j]+=x*y
    return o
def direct_score(task,mask,ncoef):
    if not strict_int(mask) or not 0<=mask<task['masks']: raise ValueError('bad candidate mask')
    bits=[0]+[int(c) for c in format(mask,f'0{task["groups"]-1}b')];a=[1.];b=[1.]
    for bit,g in zip(bits,task['roots']):
        if bit:a=conv(a,cpu_poly(g))
        else:b=conv(b,cpu_poly(g))
    if len(a)>ncoef or len(b)>ncoef:raise ValueError('candidate degree')
    product=conv(a,b)
    target=task.get('C')
    if not isinstance(target,list) or len(target)!=ncoef or any(type(x) not in (int,float) or isinstance(x,bool) or not math.isfinite(float(x)) for x in target):raise ValueError('invalid target C')
    if len(product)!=ncoef or any(not math.isfinite(x) for x in product):raise ValueError('nonfinite target convolution')
    product_error=max(abs(x-float(y)) for x,y in zip(product,target))
    a += [0.]*(ncoef-len(a));b += [0.]*(ncoef-len(b)); q=max(0.,-min(a),-min(b),max(a)-1.,max(b)-1.)
    if not math.isfinite(q):raise ValueError('nonfinite CPU convolution score')
    # Decimal root input cannot furnish an exact identity.  A larger residual
    # is an explicit numerical unknown, never a candidate acceptance.
    return {'score':q,'target_product_max_abs_error':product_error,
            'target_product_status':'numeric_consistent' if product_error<=1e-5 else 'unknown_target_mismatch'}

def validate_events(lines, expected, task, dims, input_sha, source_sha, binary_sha):
    events=[]
    for n,line in enumerate(lines,1):
        e=json.loads(line)
        if e.get('event')!='batch' or e.get('command_status')!='completed':raise ValueError(f'journal line {n}: status')
        if e.get('task')!=str(task['id']) or e.get('input_sha256')!=input_sha or e.get('source_sha256')!=source_sha or e.get('binary_sha256')!=binary_sha:raise ValueError(f'journal line {n}: binding')
        if any(e.get(k)!=v for k,v in dims.items()):raise ValueError(f'journal line {n}: dimensions')
        lo=e.get('range_start');count=e.get('range_count')
        if not strict_int(lo) or not strict_int(count) or (lo,count)!=expected[n-1]:raise ValueError(f'journal line {n}: unexpected range')
        if e.get('finite_count')!=count:raise ValueError(f'journal line {n}: nonfinite GPU result')
        inds=e.get('top_indices');scores=e.get('top_scores')
        if not isinstance(inds,list) or not isinstance(scores,list) or len(inds)!=4 or len(scores)!=4:raise ValueError(f'journal line {n}: top4 shape')
        if any(not strict_int(i) or not 0<=i<count for i in inds) or len(set(inds))!=4:raise ValueError(f'journal line {n}: indices')
        if any(type(s) not in (int,float) or isinstance(s,bool) or not math.isfinite(float(s)) for s in scores):raise ValueError(f'journal line {n}: nonfinite top score')
        candidates=[]
        for rel,gpu in zip(inds,scores):
            mask=lo+rel; cpu=direct_score(task,mask,dims['ncoef']); diff=abs(float(gpu)-cpu['score'])
            # Float32 ranking discrepancies are numerical evidence, not acceptance
            # certificates; an implausible discrepancy turns the batch into unknown.
            if diff>1e-2+1e-4*max(1.,abs(cpu['score'])):raise ValueError(f'journal line {n}: CPU/GPU mismatch {diff}')
            candidates.append({'mask':mask,'gpu_float_score':float(gpu),'cpu_direct_score':cpu['score'],'abs_difference':diff,'near_threshold_numeric_unknown':cpu['score']<1e-3,**cpu})
        events.append({'kind':'streamed_batch_numeric_ranking','coverage_credit':True,'task':task['id'],'range':[lo,lo+count],'count':count,'input_sha256':input_sha,'source_sha256':source_sha,'binary_sha256':binary_sha,'finite_count':count,'candidates_cpu_reconstructed':candidates})
    return events

def main():
 p=argparse.ArgumentParser();p.add_argument('--task',type=int,required=True);p.add_argument('--start',type=int,required=True);p.add_argument('--count',type=int,required=True);p.add_argument('--batch',type=int,choices=(65536,131072),default=131072);p.add_argument('--hours',type=float,default=4.0);p.add_argument('--out',type=Path,required=True);p.add_argument('--execute',action='store_true');p.add_argument('--smoke',action='store_true');p.add_argument('--handoff-manifest',type=Path)
 a=p.parse_args();doc=json.loads(INPUTS.read_text())
 if not 0<=a.task<len(doc['tasks']) or a.start<0 or a.count<=0 or a.hours<=0:p.error('bad task/range/cap')
 t=doc['tasks'][a.task];nf=len(doc['DC']);nc=len(doc['DC'][0]);dims={'groups':t['groups'],'nfreq':nf,'ncoef':nc}
 if t['id']!=a.task or not 1<t['groups']<32 or a.start+a.count>t['masks'] or nf>64 or nc>53:p.error('task dimensions/range invalid')
 if a.smoke and a.count % a.batch:p.error('smoke range must be an exact batch multiple')
 handoff_sha=None
 if a.execute and not a.smoke:
  if a.handoff_manifest is None:p.error('coverage execution requires --handoff-manifest')
  hand=json.loads(a.handoff_manifest.read_text());handoff_sha=sha(a.handoff_manifest)
  if hand.get('kind')!='a30-gpu52-gap-handoff-v1' or hand.get('source_inputs_sha256')!=sha(INPUTS):p.error('handoff manifest/source binding')
  if not isinstance(hand.get('layers'),list) or not hand['layers'] or len({x.get('name') for x in hand['layers'] if isinstance(x,dict)})!=len(hand['layers']):p.error('handoff requires unique nonempty layers')
  for layer in hand.get('layers',[]):
   if not isinstance(layer,dict) or sha(layer['path'])!=layer.get('sha256'):p.error('handoff layer hash mismatch')
  recomputed=a.out/'recomputed-handoff.json';cmd=[sys.executable,str(HERE/'build_handoff_manifest.py'),'--out',str(recomputed)]
  for layer in hand['layers']:cmd += ['--layer',f"{layer['name']}={layer['path']}"]
  subprocess.run(cmd,check=True,capture_output=True,text=True)
  fresh=json.loads(recomputed.read_text())
  if fresh.get('task_next')!=hand.get('task_next') or [(x['name'],x['sha256'],x.get('native_receipt_sha256')) for x in fresh['layers']]!=[(x.get('name'),x.get('sha256'),x.get('native_receipt_sha256')) for x in hand['layers']]:p.error('handoff boundaries/layers do not recompute')
  boundary=hand.get('task_next',{}).get(str(a.task))
  if type(boundary)is not int or boundary!=a.start:p.error('request must begin at exact handoff task_next')
 a.out.mkdir(parents=True,exist_ok=True);ledger=a.out/'ledger.jsonl';receipt=a.out/'receipt.json'
 if ledger.exists() or receipt.exists():p.error('refusing existing durable output')
 inp=a.out/f'task{a.task}-input.json';subprocess.run([sys.executable,str(HERE/'make_task_input.py'),'--task',str(a.task),'--out',str(inp)],check=True,capture_output=True,text=True)
 input_sha=sha(inp);source_sha=sha(SOURCE);expected=[(x,min(a.batch,a.start+a.count-x)) for x in range(a.start,a.start+a.count,a.batch)]
 manifest={'kind':'single_process_metal_stream_plan','task':a.task,'range':[a.start,a.start+a.count],'batch':a.batch,'batches':expected,'hours_cap':a.hours,'coverage_credit':not a.smoke,'input_sha256':input_sha,'source_sha256':source_sha,'handoff_manifest_sha256':handoff_sha,'dims':dims,'status':'planned' if not a.execute else 'started'};write_json(a.out/'manifest.json',manifest)
 if not a.execute:print(json.dumps(manifest,sort_keys=True));return
 remote_j=f'{REMOTE}/journal-{os.getpid()}.jsonl';remote_final=f'{REMOTE}/final-{os.getpid()}.json';local_j=a.out/'remote-stream.jsonl';local_final=a.out/'remote-final.json';status='complete';err=None;events=[];wall_started=time.monotonic();remote_source_sha=None;remote_input_sha=None;binary_sha=None
 try:
  remote(f'mkdir -p {REMOTE}')
  remote_source=f'{REMOTE}/source-{source_sha}.swift';remote_input=f'{REMOTE}/input-task{a.task}-{input_sha}.json'
  rh=remote(f"test -f {remote_source} && shasum -a 256 {remote_source} | awk '{{print $1}}' || true").stdout.strip()
  if rh!=source_sha:subprocess.run(SCP+[str(SOURCE),f'macbook-math:{remote_source}'],check=True)
  subprocess.run(SCP+[str(inp),f'macbook-math:{remote_input}'],check=True)
  remote_source_sha=remote(f"shasum -a 256 {remote_source} | awk '{{print $1}}'").stdout.strip();remote_input_sha=remote(f"shasum -a 256 {remote_input} | awk '{{print $1}}'").stdout.strip()
  if remote_source_sha!=source_sha or remote_input_sha!=input_sha:raise ValueError('strict remote byte hash mismatch')
  binary=f'{REMOTE}/native-stream-{source_sha[:16]}'
  if remote(f'test -x {binary} && echo yes || true').stdout.strip()!='yes':remote(f'TMPDIR=/Users/zhengzhihao/CodexRemoteMath/cache SWIFT_MODULECACHE_PATH=/Users/zhengzhihao/CodexRemoteMath/cache/swift swiftc -O {remote_source} -o {binary}')
  binary_sha=remote(f"shasum -a 256 {binary} | awk '{{print $1}}'").stdout.strip()
  if len(binary_sha)!=64:raise ValueError('compiled binary hash unavailable')
  ranges=';'.join(f'{x},{n}' for x,n in expected)
  # The semicolon separates ranges for Swift and must not become a remote-shell
  # command separator.  All generated range characters are decimal/commas only.
  remote(f"{binary} --input {remote_input} --ranges '{ranges}' --task {a.task} --input-sha {input_sha} --source-sha {source_sha} --binary-sha {binary_sha} --max-seconds {a.hours*3600:.6f} --journal {remote_j} --out {remote_final}",capture=False)
 except Exception as x:
  status='partial_error';err={'type':type(x).__name__,'message':str(x)}
 finally:
  try:
   subprocess.run(SCP+[f'macbook-math:{remote_j}',str(local_j)],check=True);fsync_file(local_j)
   events=validate_events(local_j.read_text().splitlines(),expected,t,dims,input_sha,source_sha,binary_sha)
   for e in events:append(ledger,e)
   if len(events)<len(expected) and status=='complete':status='partial_time_cap_or_incomplete_stream'
  except Exception as x:
   if status=='complete':status='partial_unverified_stream'
   err=err or {'type':type(x).__name__,'message':str(x)}
  try:
   subprocess.run(SCP+[f'macbook-math:{remote_final}',str(local_final)],check=True);fsync_file(local_final)
  except Exception:pass
  write_json(receipt,{'kind':'single_process_metal_stream_receipt','status':status,'task':a.task,'requested_range':[a.start,a.start+a.count],'verified_batches':len(events),'verified_prefix_end':a.start+sum(e['range'][1]-e['range'][0] for e in events),'controller_wall_seconds':time.monotonic()-wall_started,'coverage_credit':not a.smoke,'ledger_sha256':sha(ledger) if ledger.exists() else None,'remote_journal_sha256':sha(local_j) if local_j.exists() else None,'remote_source_sha256':remote_source_sha,'remote_input_sha256':remote_input_sha,'remote_binary_sha256':binary_sha,'error':err,'note':'On transport failure, remote journal is fetched and independently validated if reachable; failure to retrieve is explicitly unverified, not exclusion.'})
 print((a.out/'receipt.json').read_text())
 if status!='complete':raise SystemExit(2)
if __name__=='__main__':main()
