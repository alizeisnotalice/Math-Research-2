#!/usr/bin/env python3
"""Join slab-enumerator output to conditional slope; source must be independently audited."""
import argparse,hashlib,json,platform,sys
from fractions import Fraction as F
from pathlib import Path
from r039_endpoint_sweep import expr,slope_certificate

def convert(x):
    if isinstance(x,F): return str(x)
    if isinstance(x,dict): return {str(k):convert(v) for k,v in x.items()}
    if isinstance(x,(tuple,list)): return [convert(v) for v in x]
    return x

def main():
    p=argparse.ArgumentParser(); p.add_argument('slab'); p.add_argument('--output',required=True); p.add_argument('--bits',type=int,default=80); args=p.parse_args()
    slabpath=Path(args.slab); raw=slabpath.read_bytes(); slab=json.loads(raw)
    if slab['status']!='complete': raise ValueError('slab enumerator must report complete')
    H=F(slab['H']); delta_empty=F(slab['delta_empty'])
    if H<=0 or delta_empty<=0: raise ValueError('nonpositive slab or rejection radius')
    boxes=[tuple(expr(*e) for e in box) for box in slab['boxes']]
    if len(boxes)!=slab['counts']['active_boxes'] or len(boxes)!=len(slab['box_records']): raise ValueError('box/record mismatch')
    sweep=slope_certificate(boxes,*slab['domain'],bits=args.bits)
    delta=min(H,delta_empty,sweep['endpoint_radius'] or H)
    if delta<=0: raise ValueError('nonpositive stable radius')
    result={'scope':'conditional on completeness and validity of supplied slab enumeration; no independent enumeration audit',
      'slab_sha256':hashlib.sha256(raw).hexdigest(),'slab_code_sha256':slab['code_sha256'],
      's0':slab['s0'],'delta':str(delta),'beta_lower':str(sweep['beta_lower']),'beta_upper':str(sweep['beta_upper']),
      'tail':slab['tail'],'sweep':convert(sweep),'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
      'dependency_sha256':{name:hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() for name in ('r039_endpoint_sweep.py','r039_interval_primitives.py')},
      'python':sys.version,'platform':platform.platform(),'argv':sys.argv}
    Path(args.output).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ('s0','delta','beta_lower','beta_upper','scope')}))
if __name__=='__main__': main()
