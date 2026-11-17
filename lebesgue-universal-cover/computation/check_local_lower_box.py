"""Exact local 20-dimensional lower certificate, not a global minimum."""
from fractions import Fraction as F
from pathlib import Path
import hashlib,json,sys,runpy,io,contextlib
p=Path(__file__).with_name('013-eight-body-certificate.json');raw=p.read_bytes();r=json.loads(raw)
generator=p.with_name('certify_eight_body_upper.py')
assert hashlib.sha256(generator.read_bytes()).hexdigest()==r['script_sha256']
sys.argv=[str(generator),str(r['K']),r['target']]
with contextlib.redirect_stdout(io.StringIO()):
    rebuilt=runpy.run_path(str(generator))
assert rebuilt['H']==[tuple(q) for q in r['hull_grid_points']]
assert rebuilt['PARAMS']==r['parameter_decimals_interpreted_exactly']
K=r['K'];labels=[]
labels += [('D',j) for j in range(8*K)]
for n in (3,5,7,9):
    labels += [(f'B{n}',l,j) for l in range(n) for j in range(K+1)]
for n in (11,19):labels += [(f'P{n}',j) for j in range(n)]
labels += [('Fourier',j) for j in range(8*K)]
assert len(labels)==len(rebuilt['pts'])
lookup={q:i for i,q in enumerate(rebuilt['pts'])}
provenance=[{'point':q,'source_index':lookup[tuple(q)],'source':labels[lookup[tuple(q)]]} for q in r['hull_grid_points']]
x=list(map(F,r['parameter_decimals_interpreted_exactly']));assert len(x)==20
translation_pairs=[(0,1),(3,4),(6,7),(9,10),(12,13),(15,16),(18,19)]
assert all(x[a]**2+x[b]**2<F(1,20)**2 for a,b in translation_pairs)
assert (F(317,600)**2+F(7,60)**2)<F(3,5)**2
assert F(1,2)-F(1,80)-F(19,1200)>F(2,5)
S=r['coordinate_denominator'];V=r['hull_grid_points'];A=F(sum(a[0]*b[1]-a[1]*b[0] for a,b in zip(V,V[1:]+V[:1])),2*S*S)
assert A==F(r['rational_polygon_area'])
s=F(1,20000);eta=F(2,10**9);delta=F(13,5)*s;r0=F(7,20);alpha=(r0-delta)/(r0+eta)
assert 0<alpha<1 and (1-alpha)*r0-delta-alpha*eta==0
lower=alpha*alpha*A;assert lower>F(419,500)
print(json.dumps({'python':sys.version,'parent_certificate_sha256':hashlib.sha256(raw).hexdigest(),'generator_reexecuted':True,'hull_vertex_provenance':provenance,'box_center':[str(v) for v in x],'box_halfwidth':str(s),'dimension':20,'r0':str(r0),'eta':str(eta),'delta':str(delta),'alpha':str(alpha),'lower':str(lower),'lower_display':float(lower),'threshold':'419/500','status':'passed','scope':'one local box around the explicit F8 placement only; not all of the global domain'},indent=2))
