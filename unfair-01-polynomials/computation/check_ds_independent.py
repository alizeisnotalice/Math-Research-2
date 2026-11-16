"""Independent block audit and scattered 80-digit checks of MPS DS outputs."""
from pathlib import Path
import importlib.util,json,hashlib,os,argparse
import numpy as np
import mpmath as mp
import torch
BASE=Path(__file__).resolve().parent
parser=argparse.ArgumentParser();parser.add_argument('--source',type=Path,required=True);parser.add_argument('--roots',type=Path,required=True);parser.add_argument('--options',type=Path,required=True);parser.add_argument('--out',type=Path,required=True);args=parser.parse_args();assert Path('/Volumes/data').is_mount();SRC=args.source
spec=importlib.util.spec_from_file_location('gpu_engine',SRC);engine=importlib.util.module_from_spec(spec);spec.loader.exec_module(engine)
raw=json.loads(args.roots.read_text())['candidates'][0]
cache=json.loads(args.options.read_text());L=raw['L'];blocks=raw['blocks']
factors=[np.asarray(g['coefficients_ascending'])for g in raw['groups']]
maximum_local_relative=0.;options_checked=0
for bi,block in enumerate(cache['blocks']):
 assert (block['original_group'],block['lift_groups'])==(blocks[bi]['original_group'],blocks[bi]['lift_groups'])
 all_indices=[block['original_group']]+block['lift_groups']
 for es,options in block['options'].items():
  E=int(es);order=[]
  for option in options:
   if E<L:
    pos=option['selected_lift_positions'];assert len(set(pos))==len(pos)==E
    chosen={block['lift_groups'][q]for q in pos};order.append(sum(1<<q for q in pos))
   else:
    chosen={block['original_group']}if option['mode']=='original'else set(block['lift_groups'])
   for side,indices in [('Ablock',[j for j in all_indices if j in chosen]),('Bblock',[j for j in all_indices if j not in chosen])]:
    ref=np.array([1.])
    for j in indices:ref=np.convolve(ref,factors[j])
    got=np.asarray(option[side]['hi'],float)+np.asarray(option[side]['lo'],float)
    assert got.shape==ref.shape
    rel=float(abs(got-ref).max()/max(1.,abs(ref).max()));assert rel<1e-10
    maximum_local_relative=max(maximum_local_relative,rel)
   options_checked+=1
  if E<L:assert order==sorted(order) and len(set(order))==len(order)
  else:assert [o['mode']for o in options]==['original','all_lifts']
mp.mp.dps=80
precise=[]
for g in raw['groups']:
 key='parent_root_positive_imag'if g['type']=='original_quadratic'else'lift_root_from_positive_imag_parent'
 re,im=map(mp.mpf,g[key]);precise.append([re*re+im*im,-2*re,mp.mpf(1)])
def conv(a,b):
 c=[mp.mpf(0)]*(len(a)+len(b)-1)
 for i,x in enumerate(a):
  for j,y in enumerate(b):c[i+j]+=x*y
 return c
def independent_mask(E,code):
 mask=0
 for bi,b in enumerate(blocks):
  if E==L:opts=[1<<b['original_group'],sum(1<<j for j in b['lift_groups'])];opts=opts[:1]if bi==0 else opts
  else:
   opts=[sum(1<<b['lift_groups'][k]for k in range(L)if bits>>k&1)for bits in range(1<<L)if bits.bit_count()==E]
  mask|=opts[code%len(opts)];code//=len(opts)
 assert code==0
 return mask
rows=[];rng=np.random.default_rng(20261116)
try:
 lib=torch.mps.compile_shader(engine.METAL)
 for E in [1,4,L]:
  opts=engine.choose_options(cache,E,L);rad=[len(o)for o in opts];total=int(np.prod(rad))
  codes=np.arange(total,dtype=np.int64)if E==L else np.array(sorted(set([0,total-1]+rng.integers(0,total,16).tolist())),dtype=np.int64)
  packed=engine.prepack_options(opts);got=engine.run_batch(lib,packed,codes,rad,1e-5,materialize_coefficients=True)
  for row,code in enumerate(codes):
   mask=independent_mask(E,int(code));A=[mp.mpf(1)];B=A[:]
   for j,f in enumerate(precise):
    if mask>>j&1:A=conv(A,f)
    else:B=conv(B,f)
   refs=[]
   for a in [A,B]:
    a=[x/a[0]for x in a];a[0]=a[-1]=mp.mpf(1);refs.append(np.array([float(x)for x in a]))
   ga=engine.gpu_poly(got['Ahi'],got['Alo'],row,got['width_A'],int(got['degree_A'][row]));gb=engine.gpu_poly(got['Bhi'],got['Blo'],row,got['width_B'],int(got['degree_B'][row]))
   err=max(float(abs(a-b).max())for a,b in zip([ga,gb],refs));scale=max(1.,*(float(abs(a).max())for a in refs));assert err/scale<1e-10
   if E==L and code==0:
    H=np.asarray(raw['H_coefficients_ascending']);HL=np.zeros((len(H)-1)*L+1);HL[::L]=H
    assert max(abs(ga-H).max(),abs(gb-HL).max())<1e-10
   rows.append({'E':E,'code':int(code),'mask':mask,'max_absolute_error_vs_80digit':err,'coefficient_scale':scale,'scale_relative_error':err/scale})
 result={'verifier_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'source_sha256':hashlib.sha256(SRC.read_bytes()).hexdigest(),'options_checked':options_checked,'maximum_local_relative_error':maximum_local_relative,'random_seed':20261116,'mp_dps':80,'scattered_samples':rows,'numerical_only':True,'pid':os.getpid()}
 args.out.write_text(json.dumps(result,indent=2)+'\n')
 print(json.dumps({'options_checked':options_checked,'samples':len(rows),'max_relative':max(r['scale_relative_error']for r in rows),'pid':os.getpid()}),flush=True)
finally:
 torch.mps.synchronize();torch.mps.empty_cache()
