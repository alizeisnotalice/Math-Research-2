"""Direct Fraction AB replay of every published HGQ prefix layer, field checkpoints."""
import argparse,concurrent.futures as cf,importlib.util,json,time
from pathlib import Path
from collections import Counter
ROOT=Path(__file__).resolve().parent

def load(name):
 s=importlib.util.spec_from_file_location(name,ROOT/(name+'.py'));m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
ab=load('s234_ab');ann=load('s234_annulus')
def one(D,source,out,deadline):
 p=Path(source)/f'field-D{D}.json';target=Path(out)/f'field-D{D}.json'
 while not p.exists() and time.time()<deadline:time.sleep(2)
 if not p.exists():return {'D':D,'status':'partial_missing_input'}
 binding={q.name:ann.sha(q) for q in [Path(__file__),ROOT/'s234_ab.py']};binding['HGQ_receipt']=ann.sha(p)
 if target.exists():
  x=json.loads(target.read_text());assert x['binding']==binding and x['status']=='passed';return x
 x=json.loads(p.read_text());assert x['D']==D and x['status']=='completed_all_listed_mothers';start=time.monotonic();n=layers=0;ct=Counter()
 for row in x['records']:
  if time.time()>deadline:return {'D':D,'status':'partial_time_limit','verified_mothers':n}
  y=ab.run_core(row['row'],len(row['trace']),1000000,None)
  assert y['levels']==row['trace'],('trace mismatch',D,n,row['row'])
  assert (y['status']=='dead')==(row['status']=='dead_all_prefixes'),('terminal mismatch',D,n)
  n+=1;layers+=len(y['levels']);ct[row['status']]+=1
 assert n==x['mothers']
 z={'status':'passed','D':D,'binding':binding,'mothers':n,'layers':layers,'counts':dict(ct),'seconds':time.monotonic()-start,'scope':'Actual AB convolution at every accepted and rejected HGQ layer; no unknown treated as exclusion.'};ann.atomic(target,(json.dumps(z)+'\n').encode());return z

def main():
 if not __debug__:raise RuntimeError('no -O')
 ap=argparse.ArgumentParser();ap.add_argument('--source',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);ap.add_argument('--workers',type=int,default=2);ap.add_argument('--seconds',type=int,default=14400);args=ap.parse_args();args.out.mkdir(parents=True,exist_ok=True)
 ds=json.loads((ROOT/'row567-full/summary.json').read_text())['fields'];deadline=time.time()+args.seconds;done=[]
 with cf.ProcessPoolExecutor(max_workers=args.workers) as pool:
  jobs=[pool.submit(one,D,str(args.source),str(args.out),deadline) for D in ds]
  for job in cf.as_completed(jobs):
   z=job.result();done.append(z);print(json.dumps(z),flush=True);ann.atomic(args.out/'progress.json',(json.dumps({'records':done,'expected_fields':ds})+'\n').encode())
 z={'status':'passed' if len(done)==len(ds) and all(x['status']=='passed' for x in done) else 'partial_unknown','records':done};ann.atomic(args.out/'summary.json',(json.dumps(z)+'\n').encode())
if __name__=='__main__':main()
