#!/usr/bin/env python3
"""Only exact comparison of saved n5 clipping endpoints; no leaf reread."""
from fractions import Fraction as F
from pathlib import Path
import hashlib,json,os,sys
p=Path('computation/results/round074-borel-corners/saved-n5-n17-corners.json');d=json.loads(p.read_text());c=next(c for c in d['cases'] if c['n']==5);per_corner_budget=F((5//2)**2,25)/4;total_budget=4*per_corner_budget;lo=F(c['total_corner_deletion']['lower']);assert lo>total_budget
out={'python':sys.version,'pid':os.getpid(),'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'clipped_result_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'n':5,'per_corner_stronger_budget':str(per_corner_budget),'total_stronger_budget':str(total_budget),'total_corner_lower':str(lo),'strict_total_excess':str(lo-total_budget),'each_corner_lower_excess':{k:str(F(v['deletion']['lower'])-per_corner_budget) for k,v in c['corners'].items()},'scope':'Exact endpoint comparison only. The stronger bound each corner<=internal corner mass/4 fails for this Borel family; no finite atomic implication, no new classification or subdivision.'}
assert all(F(a)>0 for a in out['each_corner_lower_excess'].values());Path('computation/results/round074-borel-corners/stronger-borel-budget-comparison.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
