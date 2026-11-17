"""Fresh complete S1 degree8 replay: two prefix implementations and three exact root tests."""
import argparse,concurrent.futures as cf,hashlib,importlib.util,json,math,os,sys,tempfile,time
from pathlib import Path
from collections import Counter

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def digest(x):return hashlib.sha256(json.dumps(x,separators=(',',':')).encode()).hexdigest()
def module(p,name):
 s=importlib.util.spec_from_file_location(name,p);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
def atom(p,o):
 t=p.with_suffix(p.suffix+'.partial')
 with t.open('w') as f:json.dump(o,f,sort_keys=True,indent=2);f.write('\n');f.flush();os.fsync(f.fileno())
 os.replace(t,p)
def init(package):
 global PKG,P4,P7,IND,ROOT,RIND
 PKG=Path(package);F=PKG/'frozen';P4=module(F/'s1_degree8_prefix.py','p4');P7=module(F/'s1_degree8_stage5to7.py','p7');IND=module(F/'s1_degree8_independent.py','ind');ROOT=module(F/'s1_degree8_root_and_trees.py','root');RIND=module(F/'verify_s1_degree8_root_independent.py','rind')
def ascending_schur(a):
 a=list(a)
 while len(a)>1:
  g=math.gcd(*a)
  if a[-1]<0:g=-g
  a=[v//g for v in a]
  if abs(a[0])>=a[-1]:return len(a)-1
  n=len(a)-1;c,l=a[0],a[-1];a=[l*a[i]-c*a[n-i] for i in range(1,n+1)]
 return None
def independent_ascending_screen(row):
 H=RIND.H(row)
 for j,(p,q) in enumerate(RIND.POINTS):
  if sum(v*p**i*q**(16-i) for i,v in enumerate(H))<=0:return 10+j
 n=ascending_schur([v<<i for i,v in enumerate(H)])
 if n is not None:return 30+n
 n=ascending_schur([v<<i for i,v in enumerate(H[::-1])])
 return 255 if n is None else 50+n

def core_worker(expect,tmp,checkpoint,deadline):
 assert time.time()<deadline,'no success after batch deadline'
 core=tuple(expect['core']);T,N,D=core;tmp=Path(tmp);name=f'T{T}_N{N}.json';begin=time.monotonic()
 r4=P4.one(core,str(tmp/'stage4'),sha(PKG/'frozen/s1_degree8_prefix.py'),600);assert not r4.get('cached',False)
 p4=tmp/'stage4/cores'/name;raw4=json.loads(p4.read_text());assert raw4['status']=='complete'
 assert digest(raw4['records'])==expect['row4_records_sha256'] and len(raw4['records'])==expect['row4_count']
 r7=P7.one(str(p4),sha(p4),str(tmp/'stage7'),sha(PKG/'frozen/s1_degree8_stage5to7.py'),600);assert not r7.get('cached',False)
 p7=tmp/'stage7/cores'/name;raw7=json.loads(p7.read_text());assert raw7['status']=='complete'
 assert digest(raw7['records'])==expect['row7_records_sha256'] and len(raw7['records'])==expect['row7_count']
 ri=IND.stage_core(core,str(tmp/'independent'),sha(PKG/'frozen/s1_degree8_independent.py'),600);assert not ri.get('cached',False)
 pi=tmp/'independent/cores'/name;rawi=json.loads(pi.read_text());assert rawi['status']=='complete' and rawi['row4_records']==raw4['records'] and rawi['row7_records']==raw7['records']
 for key,value in r4['counts'].items():assert ri['counts'].get(key,0)==value,('stage4 counter',core,key)
 for key,value in r7['counts'].items():
  if key!='stage4_records':assert ri['counts'].get(key,0)==value,('stage7 counter',core,key)
 roots=[]
 for row in sorted({(T,N,*v[:6]) for v in raw7['records']}):
  assert ROOT.mother(row)[2]==RIND.H(row),'norm implementations disagree'
  prod=ROOT.root_screen(row);code=RIND.screen(row);asc=independent_ascending_screen(row)
  assert code==asc==prod['code'] and code!=255 and prod['status']=='rejected'
  roots.append({'row':list(row),'code':code})
 assert len(roots)==len(raw7['records']),'unexpected duplicate mother'
 result={'status':'passed','core':list(core),'row4_counts':r4['counts'],'row7_counts':r7['counts'],'row4_records':len(raw4['records']),'row7_records':len(raw7['records']),'root_records':roots,'seconds':time.monotonic()-begin};atom(Path(checkpoint)/name,result)
 for p in [p4,p7,pi]:p.unlink()
 return result

def main():
 if not __debug__:raise RuntimeError('refuse -O')
 ap=argparse.ArgumentParser();ap.add_argument('--package',type=Path,default=Path(__file__).resolve().parent);ap.add_argument('--workers',type=int,default=4);ap.add_argument('--seconds',type=int,default=14400);ap.add_argument('--out',type=Path,required=True);args=ap.parse_args();assert 1<=args.workers<=4
 if args.out.exists():raise FileExistsError('fresh output required')
 args.out.parent.mkdir(parents=True,exist_ok=True);checkpoint=args.out.with_suffix(args.out.suffix+'.checkpoints');checkpoint.mkdir();start=time.monotonic();deadline=time.time()+args.seconds;r=args.package;mh=sha(r/'core-manifest.json');vh=sha(Path(__file__))
 for x in json.loads((r/'core-manifest.json').read_text())['files']:
  p=r/x['path'];assert p.stat().st_size==x['bytes'] and sha(p)==x['sha256']
 init(str(r));data=json.loads((r/'expected.json').read_text());cores=P4.core_list();assert cores==IND.cores() and len(cores)==len(set(cores))==1140
 expected=data['cores'];assert [x['core'] for x in expected]==list(map(list,cores))
 # These exact root-filter fixtures deliberately exercise the inner branch,
 # which is not reached by any of the actual 67 rejected mothers.
 for row,code in [([-1,-1,1,0,0,0,0,0],255),([-1,-1,4,0,0,0,0,0],54)]:
  assert ROOT.root_screen(row)['code']==RIND.screen(row)==independent_ascending_screen(row)==code
 with tempfile.TemporaryDirectory(prefix='s1-degree8-full-',dir=args.out.parent) as t:
  tmp=Path(t)
  for n in ['stage4','stage7','independent']:(tmp/n/'cores').mkdir(parents=True)
  c4=Counter();c7=Counter();done=[];roots=[]
  with cf.ProcessPoolExecutor(max_workers=args.workers,initializer=init,initargs=(str(r),)) as pool:
   jobs=[pool.submit(core_worker,x,str(tmp),str(checkpoint),deadline) for x in expected]
   for job in cf.as_completed(jobs):
    x=job.result();done.append(tuple(x['core']));c4.update(x['row4_counts']);c7.update(x['row7_counts']);roots+=x['root_records']
    if len(done)%100==0 or len(done)==1140:print(json.dumps({'verified_cores':len(done),'row4':c4['row4_states'],'row7':c7['row7_states']}),flush=True)
  assert sorted(done)==cores and c4['row2_states']==3568 and c4['row3_states']==8515 and c4['row4_states']==60134
  assert c7['row5_states']==28245 and c7['row6_states']==6 and c7['row7_states']==67
  roots.sort(key=lambda x:x['row']);assert roots==data['root_records'] and len(roots)==67 and len({tuple(x['row']) for x in roots})==67
 assert sha(r/'core-manifest.json')==mh and sha(Path(__file__))==vh
 result={'status':'passed_full_S1_degree8_two_complete_prefix_implementations_and_all_root_codes','scope':'Nonrational position {1}, or reversal {7}, degree8 over any real quadratic field; arbitrary finite signed integer cofactors. No other degree8 class or general conjecture is settled.','cores':1140,'cached_cores':0,'row4_counts':dict(c4),'row7_counts':dict(c7),'unique_mothers':67,'root_codes':dict(Counter(x['code'] for x in roots)),'root_survivors':0,'actual_tree_calls':0,'root_fixture_pass_and_inner_rejection':True,'all_records_recomputed_and_compared':True,'workers':args.workers,'seconds':time.monotonic()-start,'python':sys.version,'core_manifest_sha256':mh,'verifier_sha256':vh};atom(args.out,result);print(json.dumps(result),flush=True)
if __name__=='__main__':main()
