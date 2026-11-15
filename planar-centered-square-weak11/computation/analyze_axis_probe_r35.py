#!/usr/bin/env python3
"""Enclose sampled secant slopes and identify disjoint adjacent enclosures."""
import json,sys
from fractions import Fraction as F
from pathlib import Path
from periodic_halo import sha
p=Path(sys.argv[1]);m=json.loads((p/'manifest.json').read_text());results=json.loads((p/'results.json').read_text()); assert len(results)==len(m['tasks']) and not any('error' in r for r in results)
r={x['id']:x for x in results}; groups={}
for t in m['tasks']:
 for alias in t['combinations']:
  groups.setdefault(tuple(alias['quadrants']),[]).append((F(alias['step']),r[t['id']]))
out={}
for name,pts in groups.items():
 pts.sort(key=lambda x:x[0]);segs=[]
 for (a,u),(b,v) in zip(pts,pts[1:]):
  h=b-a; low=(F(v['inner_ratio'])-F(u['outer_ratio']))/h; high=(F(v['outer_ratio'])-F(u['inner_ratio']))/h
  segs.append({'left':str(a),'right':str(b),'slope_lower':str(low),'slope_upper':str(high),'slope_lower_float':float(low),'slope_upper_float':float(high)})
 disjoint=[]
 for a,b in zip(segs,segs[1:]):
  if F(a['slope_upper'])<F(b['slope_lower']) or F(b['slope_upper'])<F(a['slope_lower']):disjoint.append({'left_interval':[a['left'],a['right']],'right_interval':[b['left'],b['right']],'left_slope':[a['slope_lower'],a['slope_upper']],'right_slope':[b['slope_lower'],b['slope_upper']]})
 out[str(name)]={'sample_points':len(pts),'segments':segs,'disjoint_adjacent_secants':disjoint}
report={'manifest_sha256':sha(p/'manifest.json'),'results_sha256':sha(p/'results.json'),'script_sha256':sha(__file__),'scope':'Enclosed secant slopes between sampled full-halo ratios, per unit normalized shift; disjoint adjacent secant enclosures imply non-affinity on their combined interval, not an exact wall location or monotonicity throughout.','directions':out}
(p/'secant_analysis.json').write_text(json.dumps(report,indent=2)+'\n')
for k,v in out.items():print(k,'points',v['sample_points'],'disjoint',len(v['disjoint_adjacent_secants']))
