#!/usr/bin/env python3
"""MPS batched numerical exploration of binary polynomial products.

Endpoints=1; a1 is kept in [.1,.9], b1=1-a1. Other coefficients are
projected to [0,1], allowing unknown support gaps. All output is numerical;
no infeasibility conclusion or exact counterexample is certified here.
"""
import argparse,json,time,os,gc,platform,hashlib
from pathlib import Path
import numpy as np
import torch

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--m',type=int,default=18);ap.add_argument('--n',type=int,default=30)
 ap.add_argument('--starts',type=int,default=16384);ap.add_argument('--steps',type=int,default=2000)
 ap.add_argument('--restart-every',type=int,default=0);ap.add_argument('--seed',type=int,default=20261115);ap.add_argument('--out',type=Path,required=True)
 a=ap.parse_args();assert Path('/Volumes/data').is_mount();a.out.mkdir(parents=True,exist_ok=True)
 assert torch.backends.mps.is_available();torch.set_num_threads(1);torch.manual_seed(a.seed)
 dev=torch.device('mps');m,n=a.m,a.n;B=a.starts
 # Independent seed rows; random positive coefficients may project exactly to 0.
 coeff=torch.nn.Parameter(torch.rand((B,m+n-4),device=dev)*.3)
 theta=torch.nn.Parameter(torch.randn((B,1),device=dev))
 one=torch.ones((B,1),device=dev)
 incidence=torch.zeros(((m+1)*(n+1),m+n+1),device=dev)
 # Build convolution incidence on CPU once, then copy to GPU.
 ci=np.zeros(((m+1)*(n+1),m+n+1),dtype=np.float32)
 for i in range(m+1):
  for j in range(n+1):ci[i*(n+1)+j,i+j]=1
 incidence=torch.tensor(ci,device=dev)
 opt=torch.optim.Adam([coeff,theta],lr=.035)
 metadata={**vars(a),'out':str(a.out),'torch':torch.__version__,'numpy':np.__version__,
           'device':'mps','platform':platform.platform(),'dtype':'float32','pid':os.getpid(),
           'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
           'scope':'degrees fixed; zeros allowed; 0.1<a1<0.9; projected numerical search only'}
 (a.out/'config.json').write_text(json.dumps(metadata,indent=2)+'\n')
 start=time.monotonic();history=[]
 def factors():
  t=.1+.8*torch.sigmoid(theta)
  aa=torch.cat((one,t,coeff[:,:m-2],one),dim=1)
  bb=torch.cat((one,1-t,coeff[:,m-2:],one),dim=1)
  return aa,bb
 try:
  for step in range(a.steps+1):
   if a.restart_every and step and step%a.restart_every==0:
    (a.out/f'restart_{step:07d}.json').write_text(json.dumps({'history':history[-1:], 'cases':cases},indent=2)+'\n')
    with torch.no_grad():coeff.uniform_(0,.3);theta.normal_()
    opt.state.clear()
   opt.zero_grad(set_to_none=True)
   aa,bb=factors()
   conv=(aa[:,:,None]*bb[:,None,:]).reshape(B,-1)@incidence
   target=(conv.detach()>=.5).to(conv.dtype)
   residual=conv-target
   per=(residual*residual).mean(dim=1)
   loss=per.sum();loss.backward();opt.step()
   with torch.no_grad():coeff.clamp_(0,1)
   if step%100==0 or step==a.steps:
    torch.mps.synchronize()
    with torch.no_grad():
     scores,ix=torch.topk(per,16,largest=False)
     aa,bb=factors();va=aa[ix].detach().cpu().numpy();vb=bb[ix].detach().cpu().numpy()
     cases=[]
     for av,bv in zip(va,vb):
      cv=np.convolve(av.astype(np.float64),bv.astype(np.float64));tar=(cv>=.5).astype(np.float64)
      cases.append({'A':av.tolist(),'B':bv.tolist(),'max_residual_float64':float(np.max(np.abs(cv-tar))),
                    'l2_residual_float64':float(np.linalg.norm(cv-tar)),'status':'numerical_only'})
     record={'step':step,'elapsed_seconds':time.monotonic()-start,'min_mse_gpu':float(scores[0].cpu()),
             'best_max_residual_float64':min(c['max_residual_float64'] for c in cases),
             'allocated_bytes':torch.mps.current_allocated_memory(),'driver_bytes':torch.mps.driver_allocated_memory()}
     history.append(record);(a.out/'checkpoint.json').write_text(json.dumps({'history':history,'cases':cases},indent=2)+'\n')
     print(json.dumps(record),flush=True)
    if step and step%500==0:
     torch.save({'coeff':coeff.detach().cpu(),'theta':theta.detach().cpu(),'step':step,
                 'optimizer':opt.state_dict()},a.out/'state.pt')
  (a.out/'result.json').write_text(json.dumps({'complete':True,'history':history,'cases':cases,
      'note':'Floating-point search is neither an exclusion certificate nor an exact counterexample.'},indent=2)+'\n')
 finally:
  # This process owns all tensors; release allocator before process exit.
  del opt,coeff,theta,one,incidence
  gc.collect();torch.mps.synchronize();torch.mps.empty_cache()
  (a.out/'release.json').write_text(json.dumps({'cleanup_called':True,'process_will_exit':True,
     'remaining_tensor_bytes_before_exit':torch.mps.current_allocated_memory()})+'\n')

if __name__=='__main__':main()
