import argparse,importlib.util,json,hashlib
from collections import Counter
from concurrent.futures import ProcessPoolExecutor,as_completed
from pathlib import Path
HERE=Path(__file__).resolve().parent
from check_union import check
def load():
 s=importlib.util.spec_from_file_location('a',HERE/'audited_chunk_verifier.py');m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
def ev(p):return [json.loads(x) for x in Path(p).read_text().splitlines() if x]
def one(root,x):
 a=load();e=[];r=a.verify_chunk(Path(root)/f"t{x['total']:02d}",x,6,e);return {'job':r['job'],'errors':e,'counts':dict(r['counts']),'verified_rows':r['verified_rows']}
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--old',type=Path,required=True);ap.add_argument('--resume',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);ap.add_argument('--cache',type=Path,required=True);ap.add_argument('--workers',type=int,default=4);z=ap.parse_args()
 m=check(z.old,z.resume,Path('/Volumes/data/CodexTemp/T2/unfair-01-polynomials/agents/support50_macbook41/plan/support50-plan.json'),Path('/Volumes/data/CodexTemp/T2/unfair-01-polynomials/agents/gpu_partial39/support50_resume41/resume-plan.json'))
 if m['errors']: raise SystemExit(m['errors'])
 jobs=[('old',x) for x in ev(z.old/'events.jsonl') if x.get('event')=='chunk' and x['total']==50]+[('resume',x) for x in ev(z.resume/'events.jsonl') if x.get('event')=='chunk' and x['total']==50]
 z.cache.mkdir(parents=True,exist_ok=True); out=[];tot=Counter(raw=0,canonical=0,excluded=0,fair=0,unknown=0); todo=[]
 binding=hashlib.sha256(json.dumps(m,sort_keys=True).encode()).hexdigest()
 for src,x in jobs:
  c=z.cache/f"{src}-{x['m']:02d}-{x['n']:02d}-{x['index']:06d}.json"
  try:
   q=json.loads(c.read_text()); assert q['binding']==binding and q['gzip']==x['gzip_sha256'] and not q['result']['errors']; out.append(q['result']);tot.update(q['result']['counts']);continue
  except Exception: todo.append((src,x))
 with ProcessPoolExecutor(max_workers=z.workers) as p:
  fs={p.submit(one,str(z.old if src=='old' else z.resume),x):(src,x) for src,x in todo}
  for f in as_completed(fs):
   src,x=fs[f];r=f.result();r['source']=src;out.append(r);tot.update(r['counts']);c=z.cache/f"{src}-{x['m']:02d}-{x['n']:02d}-{x['index']:06d}.json";t=c.with_suffix('.tmp');t.write_text(json.dumps({'binding':binding,'gzip':x['gzip_sha256'],'result':r},sort_keys=True));t.replace(c)
 out.sort(key=lambda x:(x['source'],x['job']));r={'total':50,'binding':m,'jobs':len(out),'counts':dict(tot),'errors':sum(bool(x['errors']) for x in out),'results':out};z.out.write_text(json.dumps(r,indent=2,sort_keys=True)+'\n');print({'jobs':len(out),'errors':r['errors'],'counts':r['counts']});raise SystemExit(2 if r['errors'] else 0)
if __name__=='__main__':main()
