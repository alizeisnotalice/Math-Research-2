import csv, hashlib, json
from decimal import Decimal as D, localcontext
from pathlib import Path
HERE=Path(__file__).resolve().parent
LOW=HERE/'low-batch-n11-root-summary.json'
ROOTS=[(0,(0,0)),(1,(1,0)),(12,(1,1))]; N=121; P=262144; R=256; L=16; T=12
def one(root,xy):
 p=HERE/f'n11-highp-root{root}-p{P}-r{R}-l{L}.csv'; m=json.loads(Path(str(p)+'.metadata.json').read_text())
 with p.open() as f: rows=list(csv.DictReader(f))
 assert tuple(int(m[k]) for k in ('n','N','root','P','R','L','threads'))==(11,N,root,P,R,L,T) and tuple(m['root_xy'])==xy
 assert len(rows)==R and {int(x['rep']) for x in rows}==set(range(R))
 base=int(m['seed_base']); zs=[]
 with localcontext() as c:
  c.prec=90
  for x in rows:
   q=int(x['rep']); assert (int(x['n']),int(x['N']),int(x['root']),int(x['L']))==(11,N,root,L)
   assert int(x['seed'])==(base+root*1000000007+q*1009+L*1000003)%(1<<64)
   assert int(x['proposals'])==P*L*(N-2) and int(x['dead'])==0 and 0<=int(x['accept'])<=int(x['proposals'])
   z=D(x['Z']); assert z.is_finite() and z>0; zs.append(z)
  total=sum(zs); mean=total/R; se=(sum((z-mean)**2 for z in zs)/D(R*(R-1))).sqrt()
  return {'root':root,'xy':xy,'file':p.name,'csv_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'metadata':m,'mean_h_root':str(mean),'delta_SE':str(se),'replicate_weight_ESS':str(total*total/sum(z*z for z in zs)),'max_replicate_weight_share':str(max(zs)/total),'minimum_particle_parent_ESS':min(float(x['min_parent_ess']) for x in rows),'accepted':sum(int(x['accept']) for x in rows),'proposals':sum(int(x['proposals']) for x in rows)}
m=json.loads((HERE/'run-manifest.json').read_text())
assert hashlib.sha256((HERE/'fixed_root_smc.cpp').read_bytes()).hexdigest()==m['source_sha256']
# The historical binary is intentionally excluded. Its SHA is a recorded
# value in the run manifest, not a signature and not re-verifiable here.
out={'scope':'Strict finite SMC audit; delta SE is not a certified interval or inequality proof.','source_hash_verified':True,'historical_binary_hash_recorded_only':m['binary_sha256'],'runs':[one(r,xy) for r,xy in ROOTS],'edges':[],'cross_batch_comparison':[]}
by={x['root']:x for x in out['runs']}
with localcontext() as c:
 c.prec=90
 for a,b in ((0,1),(1,12)):
  A,B=by[a],by[b]; ma,mb=D(A['mean_h_root']),D(B['mean_h_root']); sa,sb=D(A['delta_SE']),D(B['delta_SE']); ratio=mb/ma; se=ratio*((sa/ma)**2+(sb/mb)**2).sqrt()
  out['edges'].append({'from_root':a,'to_root':b,'changed_coordinate_k':0,'ratio':str(ratio),'delta_SE':str(se),'C_for_C_n_squared_over_kplus1_squared':str(ratio/D(N)),'C_delta_SE':str(se/D(N)),'strong_bound':'100','ratio_over_strong_bound':str(ratio/D(100)),'ratio_over_strong_bound_delta_SE':str(se/D(100))})
 low=json.loads(LOW.read_text()); lby={x['root']:x for x in low['runs']}
 for r in (0,1,12):
  hi,lo=by[r],lby[r]; mh,ml=D(hi['mean_h_root']),D(lo['mean_h_root']); sh,sl=D(hi['delta_SE']),D(lo['delta_SE']); q=mh/ml; qse=q*((sh/mh)**2+(sl/ml)**2).sqrt()
  out['cross_batch_comparison'].append({'root':r,'high_over_low_mean':str(q),'independent_delta_SE':str(qse),'low_file':lo['file'],'high_file':hi['file']})
(HERE/'highp-audit.json').write_text(json.dumps(out,indent=2)+'\n'); print(json.dumps(out,indent=2))
