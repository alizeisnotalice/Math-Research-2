"""Exact robust obstruction derived from the independently checked angle witnesses."""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json,sys
p=Path(__file__).parent;c=json.loads((p/'005-obstruction-certificate.json').read_text());h=json.loads((p/'004-remote-certificate.json').read_text())
assert sorted(x['cell'] for x in c['rows'])==list(range(180))
r=F(17,10000);ratios=[];slacks=[]
for row in c['rows']:
 w=list(map(F,row['weights']));L=sum(a*(abs(q[0])+abs(q[1])) for a,q in zip(w,row['normals']));g=F(row['gap_upper']);assert L>0 and g<0
 ratios.append(-g/L);slacks.append(g+r*L);assert slacks[-1]<0
radius=min(ratios)
e=F(h['epsilon'])+r
area_upper=F(h['rational_polygon_area'])+F(h['perimeter_upper'])*e+F(22,7)*e*e
assert area_upper<F(843,1000)
print(json.dumps({'python':sys.version,'source_certificate_sha256':hashlib.sha256((p/'005-obstruction-certificate.json').read_bytes()).hexdigest(),'radius':'17/10000','certified_cells':180,'witness_radius_threshold':str(radius),'threshold_display':float(radius),'largest_robust_gap':str(max(slacks)),'enlarged_cover_area_upper':str(area_upper),'area_display':float(area_upper),'scope':'C4+(17/10000)B and every subset are not universal; no global area lower bound','status':'passed'},indent=2))
