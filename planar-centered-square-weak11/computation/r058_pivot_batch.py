#!/usr/bin/env python3
"""One bounded independently reproducible rational-grid search shard."""
import argparse,hashlib,json,math,os,platform,random,subprocess,sys,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];TMP=Path('/Volumes/data/CodexTemp/T2/planar-centered-square-weak11/r058-pivot-search');BIN=TMP.parent/'r058_pivot_decomposition'

def samples(n,seed,count):
 rng=random.Random(seed);Q=10000
 for case in range(count):
  kind=case%6;alpha=(.2,.5,1,2,10)[case%5];ww=[1]*n if case%7==0 else [max(1,round(300*rng.gammavariate(alpha,1))) for _ in range(n)]
  if kind in (0,1):
   step=rng.uniform(.15,1.4)/math.sqrt(n);xy=[(j*step,(j%2)*step*rng.uniform(.1,.8) if kind==1 else 0)for j in range(n)]
  elif kind==2:
   rad=rng.uniform(.25,1.4);xy=[(rad*math.cos(2*math.pi*j/n),rad*math.sin(2*math.pi*j/n))for j in range(n)]
  elif kind==3:
   scale=rng.uniform(.8,2.5);xy=[(rng.random()*scale,rng.random()*scale)for _ in range(n)]
  elif kind==4:
   cols=max(2,round(math.sqrt(n)));step=rng.uniform(.1,.65);xy=[((j%cols)*step+rng.uniform(-.1,.1)*step,(j//cols)*step+rng.uniform(-.1,.1)*step)for j in range(n)]
  else:
   clusters=rng.choice((2,3,4));gap=rng.uniform(.2,1.2);micro=rng.uniform(.002,.08);xy=[((j%clusters)*gap+rng.uniform(-micro,micro),rng.uniform(-micro,micro))for j in range(n)]
  pts=[[round(x*Q),round(y*Q),int(w)]for (x,y),w in zip(xy,ww)]
  yield {'case':case,'kind':('chain','zigzag','ring','cloud','grid','cluster-chain')[kind],'n':n,'coordinate_denominator':Q,'points_xy_integer_weights':pts}

def main():
 a=argparse.ArgumentParser();a.add_argument('--n',type=int,required=True);a.add_argument('--seed',type=int,required=True);a.add_argument('--count',type=int,default=800);a.add_argument('--max-delete',type=int,default=3);a.add_argument('--prepare-only',action='store_true');args=a.parse_args()
 assert os.path.ismount('/Volumes/data');TMP.mkdir(parents=True,exist_ok=True);tag=f'n{args.n}-s{args.seed}-k{args.max_delete}'
 inp=TMP/(tag+'.input');out=TMP/(tag+'.output');cases=list(samples(args.n,args.seed,args.count));text=''.join(str(c['case'])+' '+str(c['n'])+' '+str(c['coordinate_denominator'])+' '+' '.join(str(v)for p in c['points_xy_integer_weights']for v in p)+'\n'for c in cases);inp.write_text(text)
 if args.prepare_only:print(inp);return
 start=time.perf_counter();r=subprocess.run([str(BIN),str(inp),str(out),str(args.max_delete)],capture_output=True,text=True);elapsed=time.perf_counter()-start
 if r.returncode:raise RuntimeError({'exit':r.returncode,'stderr':r.stderr})
 records=[]
 for line in out.read_text().splitlines():
  t=line.split();i,n=int(t[0]),int(t[1]);parts=[]
  for j in range(int(t[8])):parts.append({'mask':int(t[9+3*j]),'upper':float(t[10+3*j]),'bound_kind':int(t[11+3*j])})
  records.append({'input':cases[i],'score':float(t[2]),'deleted_mask':int(t[3]),'full_rank':int(t[4]),'eligible_subsets':int(t[5]),'actual_halo_float_if_computed':float(t[6]),'deletions_checked':int(t[7]),'leafs':parts})
 assert len(records)==args.count;records.sort(key=lambda q:-q['score']);leaders=records[:12]
 result={'scope':'Bounded search of specified <=k full-original pivot caps + eligible-component leaves, with rank and four true-prefix proxy-union fallback upper bounds. All eligibility decisions exact on rational input; rectangle areas floating numerical evidence only.','n':args.n,'seed':args.seed,'count':args.count,'max_delete':args.max_delete,'elapsed_seconds':elapsed,'exit_code':r.returncode,'score_gt_3_5_count':sum(q['score']>3.5 for q in records),'leaders':leaders,'input_sha256':hashlib.sha256(inp.read_bytes()).hexdigest(),'output_sha256':hashlib.sha256(out.read_bytes()).hexdigest(),'python':sys.version,'platform':platform.platform(),'driver_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'kernel_source_sha256':hashlib.sha256(Path(__file__).with_name('r058_pivot_decomposition.cpp').read_bytes()).hexdigest(),'kernel_binary_sha256':hashlib.sha256(BIN.read_bytes()).hexdigest(),'bound_kind':{'0':'Bmax(M,d_bbox²)','1':'eligible rank times mass','2':'mass-ascending prefix proxy union','3':'mass-descending prefix proxy union','4':'x-ascending prefix proxy union','5':'y-ascending prefix proxy union'}}
 dest=ROOT/'computation/results/r058-pivot-search';dest.mkdir(parents=True,exist_ok=True);path=dest/(tag+'.json');path.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({'path':str(path),'exit':r.returncode,'count':args.count,'elapsed':elapsed,'above3.5':result['score_gt_3_5_count'],'best':leaders[0]['score'],'best_actual':leaders[0]['actual_halo_float_if_computed'],'sha256':hashlib.sha256(path.read_bytes()).hexdigest()}),flush=True)
if __name__=='__main__':main()
