import argparse,importlib.util,json
from collections import Counter
from concurrent.futures import ProcessPoolExecutor,as_completed
from pathlib import Path
HERE=Path(__file__).resolve().parent
def load():
 s=importlib.util.spec_from_file_location('a',HERE/'audited_chunk_verifier.py');m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
def ev(p):return [json.loads(x) for x in Path(p).read_text().splitlines() if x]
def one(root,x):
 a=load();e=[];r=a.verify_chunk(Path(root)/f"t{x['total']:02d}",x,6,e);return {'job':r['job'],'errors':e,'counts':dict(r['counts']),'verified_rows':r['verified_rows']}
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--old',type=Path,required=True);ap.add_argument('--resume',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);ap.add_argument('--workers',type=int,default=4);z=ap.parse_args()
 jobs=[('old',x) for x in ev(z.old/'events.jsonl') if x.get('event')=='chunk' and x['total']==49]+[('resume',x) for x in ev(z.resume/'events.jsonl') if x.get('event')=='chunk' and x['total']==49]
 out=[];tot=Counter(raw=0,canonical=0,excluded=0,fair=0,unknown=0)
 with ProcessPoolExecutor(max_workers=z.workers) as p:
  fs={p.submit(one,str(z.old if src=='old' else z.resume),x):(src,x) for src,x in jobs}
  for f in as_completed(fs):
   src,x=fs[f];r=f.result();r['source']=src;out.append(r);tot.update(r['counts'])
 out.sort(key=lambda x:(x['source'],x['job']));r={'total':49,'jobs':len(out),'counts':dict(tot),'errors':sum(bool(x['errors']) for x in out),'results':out};z.out.write_text(json.dumps(r,indent=2,sort_keys=True)+'\n');print({'jobs':len(out),'errors':r['errors'],'counts':r['counts']});raise SystemExit(2 if r['errors'] else 0)
if __name__=='__main__':main()
