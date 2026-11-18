"""Portable independent coverage/binding audit; not another semantic replay."""
import argparse,gzip,hashlib,itertools,json,sys,time
from pathlib import Path
from collections import Counter
def sha(b):return hashlib.sha256(b).hexdigest()
def main():
 if not __debug__:raise RuntimeError('no -O')
 p=argparse.ArgumentParser();p.add_argument('--package',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();r=a.package;start=time.monotonic();assert not a.out.exists(),'fresh output required';manifest=json.loads((r/'package-manifest.json').read_text());inputs={}
 assert manifest['coverage']=={'total':54,'cap':6,'chunks':101583,'old':10784,'new':90799}
 for x in manifest['inputs']:
  q=r/x['path'];assert not q.is_symlink();b=q.read_bytes();assert len(b)==x['bytes'] and sha(b)==x['sha256'];raw=gzip.decompress(b);assert sha(raw)==x['raw_sha256'] and x['name'] not in inputs;inputs[x['name']]=raw
 read=lambda n:json.loads(inputs[n]);progress=read('progress54.json');binding=progress['binding']
 names={n:n for n in ['new-exit.json','old-exit.json','plan.json','frozen_support.py','new-events.jsonl','old-payload.py','old-events.jsonl','union-manifest.json','complement.json','old-launch.json']};names.update({'runner':'runner.py','source':'frozen_support.py','verifier':'audited_chunk_verifier.py'})
 assert set(binding)==set(names)|{'total','cap'} and binding['total']==54 and binding['cap']==6
 for k,n in names.items():assert binding[k]==sha(inputs[n])
 assert binding['source']=='d68f1010c61ad1d911749742fef8c790c0a8064bf34512b402e7e65e4703e710' and binding['verifier']=='a557bf120c6bed03adbd8b2ce153a88bf76e678ed94c6bc8815d14dc33edbe4b' and binding['runner']=='14bb62fc3f2d0bafe7661af0c10e157d3130028e16df732c5eb43bae71173283'
 expected={}
 for m in range(1,28):
  i=0
  for size in range(2,min(m+1,6)+1):
   for inside in itertools.combinations(range(1,m),size-2):expected[m,i]=(0,*inside,m);i+=1
 assert len(expected)==101583
 plan=read('plan.json');actual={}
 for j in plan['jobs']:
  if j['total']!=54:continue
  k=j['m'],j['index'];assert k not in actual and j['n']==54-j['m'];actual[k]=tuple(j['I'])
 assert actual==expected and plan['job_total']==len(plan['jobs'])==185264
 events={};orig=Counter();newkeys=set()
 for label in ('old','new'):
  for line in inputs[label+'-events.jsonl'].splitlines():
   e=json.loads(line)
   if e.get('event')!='chunk' or e.get('total')!=54:continue
   k=e['m'],e['index'];assert k in expected and k not in events and tuple(e['I'])==expected[k] and e['n']==54-k[0];events[k]=e;orig[label]+=1
   if label=='new':newkeys.add(k)
 assert events.keys()==expected.keys() and orig=={'old':10784,'new':90799}
 comp=read('complement.json');assert comp['missing_total54_jobs']==90799 and {(j['m'],j['index']):tuple(j['I']) for j in comp['missing_jobs']}=={k:expected[k] for k in newkeys}
 old=read('old-exit.json');new=read('new-exit.json');ol=read('old-launch.json')
 assert old['returncode']==255 and not old['selected_set_complete'] and old['chunks_persisted']==94465
 assert new['status']=='complete' and new['remaining']==0 and new['failures']==[] and new['persisted_chunks_total']==90799 and new['submitted']==new['finished']==90798
 assert ol['plan_sha256']==sha(inputs['plan.json']) and ol['payload_sha256']==sha(inputs['old-payload.py']) and ol['frozen_support_sha256']==binding['source']
 assert progress['status']=='complete' and progress['expected']==progress['verified']==101583 and progress['remaining']==0 and progress['failures']==[]
 end=read('exit54.json');launch=read('launch54.json');assert end['returncode']==0 and end['child_waited'] and end['remaining']==0 and end['status']=='complete' and end['verified']==101583 and end['pid']==launch['pid']==progress['pid'] and launch['runner_sha256']==binding['runner']
 seen=set();counts=Counter();archives=set()
 for x in manifest['caches']:
  assert x['path'] not in archives;archives.add(x['path']);q=r/x['path'];assert not q.is_symlink();raw=q.read_bytes();assert len(raw)==x['bytes'] and sha(raw)==x['sha256'];number=0
  for line in gzip.decompress(raw).splitlines():
   z=json.loads(line);raw=z['utf8'].encode();assert sha(raw)==z['sha256'];rel=Path(z['relativepath']);assert not rel.is_absolute() and '..' not in rel.parts and len(rel.parts)==2;k=int(rel.parent.name[1:]),int(rel.stem[1:]);assert k in expected and k not in seen;seen.add(k)
   item=json.loads(raw);e=events[k];assert item['status']=='passed' and item['errors']==[] and item['binding']==binding and item['event_sha256']==sha(json.dumps(e,sort_keys=True).encode())
   result=item['result'];assert result['counts']==e['counts'] and result['verified_rows']==e['counts']['canonical'] and result['raw_support_count']==e['counts']['raw'] and result['job']==f'm{k[0]:02}n{54-k[0]:02}/i{k[1]:06}';counts.update(result['counts']);number+=1
  assert number==x['records']
 assert seen==expected.keys() and len(archives)==27 and counts['unknown']==0
 rawpath=r.parent/'audits/raw54.json';audit=json.loads(rawpath.read_text());assert audit['status']=='passed' and audit['chunks']==101583 and audit['counts']==dict(counts) and audit['gzip_bytes']==6456579976 and audit['raw_bytes']==93647410052
 sem=json.loads((r.parent/'audits/semantic54.json').read_text());assert sem['status']=='passed' and sem['binding']==binding and sem['counts']==dict(counts) and sem['raw_audit_sha256']==sha(rawpath.read_bytes()) and sem['remote_exit_sha256']==sha(inputs['exit54.json'])
 assert json.loads((r.parent/'audits/semantic54-exit.json').read_text())=={'returncode':0,'child_waited':True}
 result={'status':'passed_full_archive_coverage_and_binding','chunks':len(seen),'counts':dict(counts),'origin':dict(orig),'raw_audit_sha256':sha(rawpath.read_bytes()),'manifest_sha256':sha((r/'package-manifest.json').read_bytes()),'auditor_sha256':sha(Path(__file__).read_bytes()),'seconds':time.monotonic()-start,'python':sys.version,'scope':'Exact complete total54 receipt coverage and original-byte audit bindings; archive audit is not a second semantic replay.'}
 a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
if __name__=='__main__':main()
