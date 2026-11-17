"""Checkpointed per-field S234 HGQ necessary trees; unknowns retained explicitly."""
import argparse,concurrent.futures as cf,importlib.util,json,time
from pathlib import Path
from collections import Counter
ROOT=Path(__file__).resolve().parent
s=importlib.util.spec_from_file_location('hgq',ROOT/'s234_hgq.py');h=importlib.util.module_from_spec(s);s.loader.exec_module(h);a=h.a

def one(D,screen,out,deadline):
 source=Path(screen)/f'field-D{D}.json';target=Path(out)/f'field-D{D}.json'
 while not source.exists() and time.time()<deadline:time.sleep(2)
 if not source.exists():return {'D':D,'status':'partial_missing_input'}
 binding={p.name:a.sha(p) for p in [Path(__file__),ROOT/'s234_hgq.py',ROOT/'s234_annulus.py',a.RING]};binding['screen_sha256']=a.sha(source)
 if target.exists():
  z=json.loads(target.read_text());assert z['status']=='completed_all_listed_mothers' and z['binding']==binding
  return {'D':D,'status':'complete_cached','counts':z['counts']}
 x=json.loads(source.read_text());assert x['status']=='complete' and x['D']==D
 records=[];counts=Counter();start=time.monotonic()
 for row in x['survivors']:
  if time.time()>deadline:return {'D':D,'status':'partial_time_limit','mothers_processed':len(records)}
  z=h.tree(row);records.append(z);counts[z['status']]+=1
 assert len({tuple(z['row']) for z in records})==len(records)
 result={'status':'completed_all_listed_mothers','D':D,'binding':binding,'counts':dict(counts),'records':records,'mothers':len(records),'seconds':time.monotonic()-start}
 a.atomic(target,(json.dumps(result,separators=(',',':'))+'\n').encode());return {'D':D,'status':'complete','counts':dict(counts),'mothers':len(records),'seconds':result['seconds']}

def main():
 if not __debug__:raise RuntimeError('no -O')
 ap=argparse.ArgumentParser();ap.add_argument('--screen',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);ap.add_argument('--workers',type=int,default=4);ap.add_argument('--seconds',type=int,default=14400);args=ap.parse_args();args.out.mkdir(parents=True,exist_ok=True)
 fields=json.loads((ROOT/'row567-full/summary.json').read_text())['fields'];deadline=time.time()+args.seconds;done=[];counts=Counter();start=time.monotonic()
 with cf.ProcessPoolExecutor(max_workers=args.workers) as pool:
  jobs=[pool.submit(one,D,str(args.screen),str(args.out),deadline) for D in fields]
  for job in cf.as_completed(jobs):
   z=job.result();done.append(z);counts.update(z.get('counts',{}));print(json.dumps(z),flush=True);a.atomic(args.out/'progress.json',(json.dumps({'records':done,'expected_fields':fields,'counts':dict(counts)})+'\n').encode())
 status='complete' if len(done)==len(fields) and all(z['status'].startswith('complete') for z in done) else 'partial_unknown'
 result={'status':status,'records':done,'counts':dict(counts),'seconds':time.monotonic()-start,'scope':'Only dead trees exclude; depth/state/time unknowns are not excluded.'};a.atomic(args.out/'summary.json',(json.dumps(result)+'\n').encode());print(json.dumps(result),flush=True)
if __name__=='__main__':main()
