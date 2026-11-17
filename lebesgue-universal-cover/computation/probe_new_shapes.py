"""Numerical placement screening only: finite-direction outer relaxation."""
import json,math,sys
from pathlib import Path
import numpy as np
import scipy
from scipy.optimize import linprog
from fractions import Fraction
r=json.loads((Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).with_name('004-remote-certificate.json')).read_text());V=np.array(r['hull_grid_points'],float)/r['coordinate_denominator'];e=np.roll(V,-1,axis=0)-V
N=np.column_stack((e[:,1],-e[:,0]));N/=np.linalg.norm(N,axis=1)[:,None];b=(N*V).sum(axis=1)+float(Fraction(r['epsilon']))
rows=[]
for n in (3,5,7,9,11,13,15,17,19,21):
 R=1/(2*math.cos(math.pi/(2*n)));best=None
 for j in range(180):
  theta=2*math.pi*j/(n*180);a=theta+np.arange(n)*2*math.pi/n
  P=R*np.column_stack((np.cos(a),np.sin(a)))
  h=(N@P.T).max(axis=1)
  z=linprog([0,0,1],A_ub=np.column_stack((N,-np.ones(len(N)))),b_ub=b-h,bounds=[(None,None)]*3,method='highs')
  assert z.success
  if best is None or z.fun<best['slack']:best={'theta':theta,'slack':z.fun,'translation':z.x[:2].tolist()}
 rows.append({'n':n,'best':best,'grid_halfstep_lipschitz':R*math.pi/(n*180)})
print(json.dumps({'numpy':np.__version__,'scipy':scipy.__version__,'python':sys.version,'angles':180,'results':rows,'scope':'floating LP; only vertices of regular Reuleaux bodies; candidate exterior polygon relaxation, not certificate'},indent=2))
