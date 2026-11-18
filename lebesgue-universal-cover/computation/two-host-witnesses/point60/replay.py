"""Replay finite witnesses and strict all-pose escape/area bounds."""
from pathlib import Path
from fractions import Fraction as F
import sys,json,math,subprocess,hashlib,shutil
D=Path(__file__).resolve().parent
P=D.parent
C=Path(sys.argv[1]).resolve()
subprocess.run([sys.executable,str(D/'points60.py'),str(C)],check=True)
# Run the all-angle verifier when requested and validate all eight endpoint bounds.
if '--reuse-checked' not in sys.argv[2:]:
    subprocess.run([sys.executable,'-I','-B',str(P/'verify_cert.py'),str(C)],check=True)
sys.path.insert(0,str(P))
from cross_common import q,arb,pi,axis,Host,A,ctx

# Recheck host radius, using 64 directions with rational coordinates.
seed=[]
for i in range(64):
    a=2*math.pi*i/64
    sign=1 if math.cos(a)>=0 else -1
    b=(a-(0 if sign==1 else math.pi)+math.pi)%(2*math.pi)-math.pi
    seed.append((sign,round(10**12*math.tan(b/2))))
u=[axis(s)[0] for s in seed]
assert all(x[0]*y[1]-x[1]*y[0]>0 for x,y in zip(u,u[1:]+u[:1]))
assert all(q(x[0]*y[0]+x[1]*y[1])>(2*pi/64+A('.0001')).cos() for x,y in zip(u,u[1:]+u[:1]))
cap=A('.59')*(pi/64+A('.00005')).cos()
pose155=json.loads((C/'free-p30-reoptimized-pose.json').read_text())
pose175=json.loads((C/'seven-point-families-pose.json').read_text())['pose']
for h in (Host(pose155,False),Host(pose175,True)):
    assert all(h.row(s)[2]<cap for s in seed)

# ODE comparison h''+h <= 1 gives this support sampling loss. Shrink and
# 10^-12 rounding add at most .58*10^-9 + 1.5*10^-12.
loss=1-(pi/60).cos()+A('.58')*q(F(1,10**9))+q(F(3,2*10**12))
assert loss<A('.001371001')
base={'r155-r151':F('.00154117'),'r155-r216':F('.00149681'),
      'r175-r151':F('.00145235'),'r175-r216':F('.00141054')}
strict=json.loads((P/'strict-result.json').read_bytes())
assert strict['labels_sha256']==hashlib.sha256((P/'certificate-labels.json').read_bytes()).hexdigest()
assert set(strict['lower_bounds'])=={f'{name}-{sg}' for name in base for sg in (1,-1)}
for name,threshold in base.items():
    assert all(arb(strict['lower_bounds'][f'{name}-{sg}'])>q(threshold) for sg in (1,-1))
rows={}
for name,e in base.items():
    d=e-F('.001371001')
    assert d>0 and (F(59,100)+d)**2-F(1,4)<F(8,25)**2
    area=F(25,16)*d*d
    rows[name]={'escape_lower':str(d),'area_increase_lower':str(area)}
    print(name,'escape >',d,'area increase >',area,flush=True)
obj={'sampling_loss_upper':str(loss.upper()),'rows':rows,
     'points_sha256':hashlib.sha256((D/'points60.json').read_bytes()).hexdigest(),
     'precision_bits':ctx.prec}
(D/'finite-result.json').write_text(json.dumps(obj,indent=2)+'\n')
