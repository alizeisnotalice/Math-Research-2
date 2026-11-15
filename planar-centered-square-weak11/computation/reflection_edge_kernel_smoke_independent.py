#!/usr/bin/env python3
"""Independent bounded source-window and sparse exact-area tests of new edge/corner kernels."""
from fractions import Fraction as F
from pathlib import Path
import argparse,importlib.util,hashlib,json,sys,time,os
assert os.path.ismount('/Volumes/data');sys.dont_write_bytecode=True;t0=time.perf_counter()
here=Path(__file__).resolve();ap=argparse.ArgumentParser();ap.add_argument('--kernel-dir',type=Path,default=here.parent);args=ap.parse_args();base=args.kernel_dir
sys.path.insert(0,str(base));modules={};sources={}
for key,name in [('edge','reflection_edge_states.py'),('corner','reflection_edge_corners.py')]:
    path=base/name;blob=path.read_bytes();sources[name]=dict(sha256=hashlib.sha256(blob).hexdigest(),bytes=len(blob))
    spec=importlib.util.spec_from_file_location('independent_'+key,path);module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);modules[key]=module
union=base/'periodic_halo.py';sources[union.name]=dict(sha256=hashlib.sha256(union.read_bytes()).hexdigest(),bytes=union.stat().st_size)
pts=[(F(0),F(0),F(1)),(F(1,2),F(3,4),F(2)),(F(7,4),F(11,4),F(3,10))]
win=(F(-1),F(5),F(-2),F(4));R=F(2);C=F(3);px=F(2);py=F(3)
edge_states=['bulk','left','right','bottom','top','seam-old','seam-new']
corner_states=['corner-lb','corner-lt','corner-rb','corner-rt','seam-bottom-old','seam-bottom-new','seam-top-old','seam-top-new']
def reference(state):
    out={}
    for mirror in ([False,True] if state=='seam-new' or state.endswith('-new') else [False]):
      for a,b,m in pts:
       for i in range(-12,13):
        if (state=='left' or state.startswith('corner-l')) and i<0:continue
        if (state=='right' or state.startswith('corner-r')) and i>-1:continue
        if state.startswith('seam-') and i>0:continue
        x=a+i*px
        if mirror:x=2*C-x
        for j in range(-12,13):
         if (state=='bottom' or (state.startswith('corner-') and state.endswith('b')) or state.startswith('seam-bottom-')) and j<0:continue
         if (state=='top' or (state.startswith('corner-') and state.endswith('t')) or state.startswith('seam-top-')) and j>-1:continue
         y=b+j*py
         if win[0]-R<=x<=win[1]+R and win[2]-R<=y<=win[3]+R:out[x,y]=out.get((x,y),F(0))+m
    return out
filters=[]
for key,states in [('edge',edge_states),('corner',corner_states)]:
 for state in states:
    got=modules[key].source_atoms(pts,px,py,win,R,state,C);want=reference(state);assert got==want
    filters.append(dict(state=state,atoms=len(got),mass=str(sum(got.values(),F(0)))))
# Sparse unit source grid: reflected tail at r=1 is 7/8<1, distinct centers have Linfty distance>=2.
# Thus only singleton halos exist; the following clipped areas are known exactly without bbox code.
data=dict(Px='4',Py='5',lam='1',pts=[['0','0','1']]);seam=dict(input_sha256='independent_sparse',c_edge_upper='1',rational_reflection_axis='1')
assert F(1,20)*(1+4)*(1+F(5,2))==F(7,8)<1
expected={'bulk':F(1),'left':F(5,2),'right':F(3,2),'bottom':F(5,2),'top':F(3,2),'seam-old':F(4),'seam-new':F(4),'corner-lb':F(25,4),'corner-lt':F(15,4),'corner-rb':F(15,4),'corner-rt':F(9,4),'seam-bottom-old':F(14),'seam-bottom-new':F(14),'seam-top-old':F(10),'seam-top-new':F(10)}
area_rows=[]
for key,states in [('edge',edge_states),('corner',corner_states)]:
 modules[key].INPUT_HASH='independent_sparse'
 for state in states:
    result=modules[key].calculate(data,state,seam if state.startswith('seam-') else None,bits=10)
    lo,hi=F(result['inner_area']),F(result['outer_area']);assert lo<=expected[state]<=hi
    area_rows.append(dict(state=state,exact_area=str(expected[state]),inner=str(lo),outer=str(hi),boxes=result['counts']['boxes']))
for name,record in sources.items():assert hashlib.sha256((base/name).read_bytes()).hexdigest()==record['sha256']
here=Path(__file__).resolve();result=dict(status='PASS',python=sys.version,arithmetic='Fraction; exact toy profiles plus certified kernel inner/outer intervals',random_seed=None,floating_tolerance=None,source_snapshots=sources,source_filter_tests=filters,independent_sparse_profile_tests=area_rows,script_sha256=hashlib.sha256(here.read_bytes()).hexdigest(),elapsed_seconds=time.perf_counter()-t0,scope='Small deterministic source-filter and exact sparse-area tests only; no re-execution of expensive R35/R41 profiles and no claims that their leading/corner results are complete')
path=here.parent/'results/reflection-edge-kernel-independent/smoke.json';path.parent.mkdir(parents=True,exist_ok=True);path.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
print(json.dumps({k:result[k] for k in ['status','source_snapshots','elapsed_seconds']},ensure_ascii=False,indent=2))
