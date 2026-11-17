"""Full S234 replay: complete integral box, exact prefixes, every root byte, actual AB/HGQ trees."""
import argparse,concurrent.futures as cf,gzip,hashlib,importlib.util,json,math,os,shutil,sys,tempfile,time
from collections import Counter,defaultdict
from pathlib import Path

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def module(p,n):
 s=importlib.util.spec_from_file_location(n,p);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
def atomic(p,o):
 t=p.with_suffix(p.suffix+'.partial')
 with t.open('w') as f:json.dump(o,f,sort_keys=True,indent=2);f.write('\n');f.flush();os.fsync(f.fileno())
 os.replace(t,p)
def init(package):
 global PACKAGE,R4,R7,RC,TREE,AB
 PACKAGE=Path(package);f=PACKAGE/'frozen';R4=module(f/'s234_row4.py','r4full');R7=module(f/'s234_row567.py','r7full');RC=module(f/'verify_s234_annulus_descending.py','rcfull');TREE=module(f/'s234_hgq.py','treefull');AB=module(f/'s234_ab.py','abfull')
def field(D,records,tmp,deadline,checkpoint):
 start=time.monotonic();tmp=Path(tmp);c4,rows=R4.one(D);p4=tmp/'row4'/f'field-D{D}.json';R4.atomic(p4,{'D':D,'status':'complete','counts':c4,'rows':rows});del rows
 p7=tmp/'row7'/p4.name;r7=R7.one(D,str(p4),str(p7),sha(PACKAGE/'frozen/s234_row567.py'));assert r7['D']==D;c7=r7['counts'];p4.unlink()
 rc=RC.one(D,str(tmp/'row7'),str(PACKAGE/'certificate/input/annulus'),str(tmp/'root-audits'),deadline,sha(PACKAGE/'frozen/verify_s234_annulus_descending.py'));assert rc['status']=='passed'
 root=json.loads((PACKAGE/'certificate/input/annulus'/p4.name).read_text());assert root['status']=='complete' and root['D']==D;assert [x['row'] for x in records]==root['survivors'];del root
 ct=Counter();layers=deadmax=maxlive=0
 for record in records:
  assert time.time()<deadline,'batch time limit: no success receipt'
  row=record['row'];k=record['death_layer'];assert record['status']=='dead_all_prefixes' and 1<=k<=71
  h=TREE.tree(row,depth=k,cap=25000,seconds=60);assert h['status']=='dead_all_prefixes' and h['depth_reached']==k
  ab=AB.run_core(row,depth=k,cap=25000,seconds=60);assert ab['status']=='dead' and ab['depth_reached']==k and ab['levels']==h['trace']
  digest=hashlib.sha256(json.dumps(h['trace'],sort_keys=True,separators=(',',':')).encode()).hexdigest();assert digest==record['trace_sha256']
  ml=max(x['live'] for x in h['trace']);assert ml==record['max_live'];assert h['trace'][-1]['state_sha256']==record['final_state_sha256']
  layers+=k;deadmax=max(deadmax,k);maxlive=max(maxlive,ml);ct[h['status']]+=1
 p7.unlink();R4.row4.cache_clear()
 result={'status':'passed','D':D,'row4_counts':c4,'prefix_counts':c7,'root_counts':rc['counts'],'tree_counts':dict(ct),'tree_layers':layers,'deadmax':deadmax,'maxlive':maxlive,'seconds':time.monotonic()-start,'scope':'Actual recomputation, including every descending root byte and every AB/HGQ layer.'};atomic(Path(checkpoint)/f'field-D{D}.json',result);return result

def main():
 if not __debug__:raise RuntimeError('no -O')
 ap=argparse.ArgumentParser();ap.add_argument('--package',type=Path,default=Path(__file__).resolve().parent);ap.add_argument('--workers',type=int,default=4);ap.add_argument('--seconds',type=int,default=14400);ap.add_argument('--out',type=Path,required=True);args=ap.parse_args();assert 1<=args.workers<=12
 if args.out.exists():raise FileExistsError('Refusing existing output: choose a fresh path')
 r=args.package;args.out.parent.mkdir(parents=True,exist_ok=True);checkpoint=args.out.with_suffix(args.out.suffix+'.checkpoints');checkpoint.mkdir();start=time.monotonic();deadline=time.time()+args.seconds;corehash=sha(r/'core-manifest.json');driverhash=sha(Path(__file__))
 for item in json.loads((r/'core-manifest.json').read_text())['files']:
  p=r/item['path'];assert p.stat().st_size==item['bytes'] and sha(p)==item['sha256']
 cert=r/'certificate';meta=json.loads((cert/'package-manifest.json').read_text());fields=[D for D in range(2,190) if (D<=47 or D%4==1) and all(D%(i*i) for i in range(2,math.isqrt(D)+1))];assert len(fields)==58 and meta['counts']['fields']==58
 records=defaultdict(list)
 for item in meta['tree_index_parts']:
  p=cert/item['path'];assert sha(p)==item['sha256'];data=gzip.decompress(p.read_bytes());assert hashlib.sha256(data).hexdigest()==item['raw_sha256'];number=0
  for line in data.splitlines():x=json.loads(line);assert x['D']==x['row'][0] and x['D'] in fields;records[x['D']].append(x);number+=1
  assert number==item['records']
 assert sum(map(len,records.values()))==68593 and len({tuple(x['row']) for seq in records.values() for x in seq})==68593
 with tempfile.TemporaryDirectory(prefix='s234-full-',dir=args.out.parent) as tmp:
  tmp=Path(tmp)
  for name in ['row4','row7','root-audits']:(tmp/name).mkdir()
  c4=Counter();c7=Counter();roots=Counter();trees=Counter();layers=deadmax=maxlive=0;done=[]
  with cf.ProcessPoolExecutor(max_workers=args.workers,initializer=init,initargs=(str(r),)) as pool:
   jobs=[pool.submit(field,D,records[D],str(tmp),deadline,str(checkpoint)) for D in fields]
   for job in cf.as_completed(jobs):
    z=job.result();done.append(z['D']);c4.update(z['row4_counts']);c7.update(z['prefix_counts']);roots.update(z['root_counts']);trees.update(z['tree_counts']);layers+=z['tree_layers'];deadmax=max(deadmax,z['deadmax']);maxlive=max(maxlive,z['maxlive']);print(json.dumps({'D':z['D'],'verified_fields':len(done),'tree_mothers':trees['dead_all_prefixes']}),flush=True)
  assert sorted(done)==fields and c4['ak_pass']==61596 and c4['core_pass']==5181392 and c4['row4_pass']==4459150
  assert c7['row7_states']==5827896 and c7['state_f_combinations']==66621398
  assert roots['unique_bases']==4796098 and roots['unique_mothers']==54726220 and roots['survived']==68593
  assert trees=={'dead_all_prefixes':68593} and layers==621704 and deadmax==71 and maxlive==7
 assert corehash==sha(r/'core-manifest.json') and driverhash==sha(Path(__file__)),'frozen input changed during replay'
 result={'status':'passed_full_S234_integral_box_prefix_every_root_byte_and_actual_AB_HGQ','scope':'Degree7 real quadratic P with nonrational positions {2,3,4}; arbitrary finite signed integer cofactors, including constants. General conjecture remains open.','fields':58,'row4_counts':dict(c4),'prefix_counts':dict(c7),'root_counts':dict(roots),'tree_counts':dict(trees),'tree_layers':layers,'deadmax':deadmax,'maxlive':maxlive,'direct_AB_replayed':True,'all_fields_recomputed':True,'cached_fields':0,'workers':args.workers,'seconds':time.monotonic()-start,'python':sys.version,'core_manifest_sha256':corehash,'verifier_sha256':driverhash};atomic(args.out,result);print(json.dumps(result),flush=True)
if __name__=='__main__':main()
