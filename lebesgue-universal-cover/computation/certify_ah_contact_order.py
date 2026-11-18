"""Fixed-angle constants for monotone contact graph and reflected corner exclusion."""
from pathlib import Path
from sigma_parameter import SIGMA_OUT, radians
from fractions import Fraction as F
import runpy,io,contextlib,sys,json,hashlib
from rational_trig import PI,trig,scale
base=Path(__file__).parent;argv=sys.argv;sys.argv=[str(base/'certify_ah_area.py'),'64']
with contextlib.redirect_stdout(io.StringIO()):a=runpy.run_path(str(base/'certify_ah_area.py'))
sys.argv=argv;d=a['d'];I=d['I'];V=d['V'];add=d['add'];sub=d['sub'];mul=d['mul'];dot=d['dot'];norm2=d['norm2']
D=sub(d['T'],d['Z']);assert D[0].lo>0 and D[1].lo>0
assert dot(d['dL'],D).lo>0
assert all(x.lo>0 for x in sub(d['box'],d['Nbox']))
nN=V(F(-1,2),-d['r3']/2);eN=d['dN'];T=d['T'];dist=dot(nN,T)-F(1,2)
Nt=sub(sub(T,mul(dist,nN)),mul((1-dist*dist).sqrt(),eN))
tT=dot(sub(Nt,d['M']),eN);tau=norm2(sub(d['E'],d['M'])).sqrt()
assert tT.lo>0 and tT.hi<tau.lo
Tp=a['Tp'];assert Tp[1].lo>F(518,1000)
for P in (d['Y'],d['Z'],d['T'],a['Xas']):assert P[1].lo>F(518,1000)
def unit(q, coefficient=0):
 z=radians(q, coefficient);return V(I(*trig(z,'cos')),I(*trig(z,'sin')))
def meet(u,v):
 det=u[0]*v[1]-u[1]*v[0]
 return V((v[1]-u[1])/(2*det),(u[0]-v[0])/(2*det))
np4=unit(F(3,2), -1)
D2p=meet(unit(F(4,3)),np4);D3p=meet(unit(F(5,3)),np4);D1=V(0,-I(1)/d['r3'])
assert all(P[1].hi<F(-496,1000) for P in (D1,D2p,D3p))
assert F(518,1000)-F(-496,1000)>1
# The C' cut increases strictly along the DC direction.
np5=unit(F(11,6), -1)
assert dot(np5,d['dL']).lo>0
print(json.dumps({'python':sys.version,'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),"sigma_degrees_interval": SIGMA_OUT,'t_T_interval':tT.out(),'tau_interval':tau.out(),'line_D_components':[x.out() for x in D],'line_eL_dot_D_lower':str(dot(d['dL'],D).lo),'AH_lower_y':'259/500','Dprime_upper_y':'-62/125','vertical_separation_lower':'507/500','Dprime_vertices':[[x.out() for x in P] for P in (D1,D2p,D3p)],'Cprime_cut_derivative_lower':str(dot(np5,d['dL']).lo),'passed':True,'scope':'contact-order and corner-separation constants; uses prior boundary proofs and does not by itself prove the full universal cover'},indent=2))
