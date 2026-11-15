#!/usr/bin/env python3
import argparse,hashlib,importlib.util,json,sys,platform
from collections import Counter
from fractions import Fraction as F
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--batch-dir',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
base=Path(__file__).parent;spec=importlib.util.spec_from_file_location('core',base/'finite_source_cut_certificate.py');core=importlib.util.module_from_spec(spec);spec.loader.exec_module(core)
rows=[];identities=Counter();batch_sources=[]
for j in range(12):
 path=a.batch_dir/f'batch-{j}.json';data=json.loads(path.read_text());batch_sources.append({'path':str(path),'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'bytes':path.stat().st_size})
 assert data['status']=='PASS' and data['script_sha256']=='e1d356ac0771b36e3a2145775e98ee27d3211c56447f29fca6939e2b540e9b41'
 for k,r in enumerate(data['cases']):
  inp=r['input'];pts=inp['points_quarters'];ms=inp['mass_sixteenths'];identity=json.dumps(sorted([list(p)+[m] for p,m in zip(pts,ms)]),separators=(',',':'));identities[identity]+=1
  s=r['geometric_critical_mask'];ids=[i for i in range(len(pts)) if s>>i&1];span=max(max(pts[i][d] for i in ids)-min(pts[i][d] for i in ids) for d in(0,1));mass=sum(ms[i] for i in ids)
  inn,out,_=core.rectangles(pts,ms,r['bits']);G=sum((b-a)*(d-c) for a,b,c,d,_ in out)-sum((b-a)*(d-c) for a,b,c,d,_ in inn)
  assert G>=0
  rows.append({'batch':j,'case':k,'n':len(pts),'bound':r['geometric_bound'],'raw_bound':r['raw_bound'],'critical_mask':s,'critical_rank':len(ids),'critical_span2_over_mass':str(F(span*span,mass)),'critical_compact':span*span<=mass,'proper_critical':s!=(1<<len(pts))-1,'compact_cuts':r['compact_cuts'],'noncompact_cuts':r['noncompact_cuts'],'terminal_improved_cuts':r['terminal_improved_cuts'],'label_uncertainty_area':str(F(G,1<<(2*r['bits']))),'transport_cost_error_bound':str(F(16*G,(1<<(2*r['bits']))*min(ms)))})
best=max(rows,key=lambda r:F(r['bound']));worsterror=max(rows,key=lambda r:F(r['transport_cost_error_bound']))
summary={'status':'PASS','python':sys.version,'platform':platform.platform(),'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'cases_evaluated':len(rows),'distinct_inputs_up_to_source_reordering':len(identities),'repeated_evaluations':len(rows)-len(identities),'n_counts':dict(sorted(Counter(r['n'] for r in rows).items())),'total_all_cuts':sum((1<<r['n'])-1 for r in rows),'compact_cuts':sum(r['compact_cuts'] for r in rows),'noncompact_cuts':sum(r['noncompact_cuts'] for r in rows),'terminal_improved_cuts':sum(r['terminal_improved_cuts'] for r in rows),'proper_critical_cases':sum(r['proper_critical'] for r in rows),'noncompact_critical_cases':sum(not r['critical_compact'] for r in rows),'maximum_bound_case':best,'maximum_certified_precision_error_case':worsterror,'all_bounds_below_3_3':all(F(r['bound'])<F(33,10) for r in rows),'batch_sources':batch_sources,'rows':rows,'scope':'Finite inputs at threshold1 only; exact rational precision errors and counts; no general C2 bound. Repeated controls are explicitly counted.'}
a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps({k:v for k,v in summary.items() if k not in ('rows','batch_sources')},indent=2))
