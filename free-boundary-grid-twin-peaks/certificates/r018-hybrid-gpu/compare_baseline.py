import json,subprocess,os,hashlib
from pathlib import Path
archive=Path(__file__).resolve().parent
root=Path('/Volumes/data/CodexTemp/T2/a022-r018-baseline-reproduction');root.mkdir(parents=True,exist_ok=True)
src=archive.parent/'r017-bitset/rank_audit_final_source.swift'
new=archive
env=os.environ.copy();env.update(TMPDIR=str(root),TMP=str(root),TEMP=str(root),CLANG_MODULE_CACHE_PATH='/Volumes/data/CodexTemp/cache/swift-modules')
subprocess.run(['swiftc','-O','-module-cache-path',env['CLANG_MODULE_CACHE_PATH'],str(src),'-o',str(root/'old_rank')],check=True,env=env)
results=[]
for tau in range(2,11):
 f=new/f'strat-n3-t{tau}-l0.json'
 data=json.loads(f.read_text()); assert data['n']==3 and data['L']==0
 seed=data['seed_hex'];P=data['P'];oldfile=root/f'old-n3-t{tau}.jsonl'
 subprocess.run([str(root/'old_rank'),'3',str(tau),str(P),seed,'strat',str(oldfile)],check=True,env=env,stdout=subprocess.DEVNULL)
 rows=[x for line in oldfile.read_text().splitlines() if (x:=json.loads(line)).get('kind')=='rank']
 states=[[] for _ in range(P)];checks=0
 for x in data['traces']:
  r=x['rank']; y=rows[r-1]
  assert y['rank']==r and x['allowed']==y['allowed']
  assert x['cdf'][-1]==y['total']==y['cdf'][-1]
  if x.get('dead'):
   assert y['dead'];break
  if not x.get('terminal'):
   assert x['parents']==y['parents'];states=x['states']
  for i,w in enumerate(states):
   o=[0]*32;b=[0]*32;seen=set(w)
   for v in w:o[v//32]|=1<<(v%32)
   for v in seen:
    for z in ([v-1] if v%3 else [])+([v+1] if v%3<2 else [])+([v-3] if v//3 else [])+([v+3] if v//3<2 else []):
     if z not in seen:b[z//32]|=1<<(z%32)
   rr=[w[0] if w else 0,w[tau-1] if len(w)>=tau else 0]
   assert y['state_words'][66*i:66*(i+1)]==o+b+rr,(tau,r,i)
   checks+=1
 results.append({'n':3,'tau':tau,'P':P,'seed_hex':seed,'old_actual_GPU_vs_new_L0_states_checked':checks,'new_trace_sha256':hashlib.sha256(f.read_bytes()).hexdigest(),'old_trace_sha256':hashlib.sha256(oldfile.read_bytes()).hexdigest()})
out={'status':'passed','scope':'Actual frozen r017 GPU rank-audit executable vs r018 hybrid host L0; all n3 strata, stratified parents. Old rank N stores preterminal state, compared to new rank N-1 state. No cross-kernel L>0 identity claimed.','rows':results}
(root/'direct-baseline-certificate.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
