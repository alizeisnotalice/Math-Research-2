"""Replay the scale-interval transfer; the .83745 native area trees are inherited inputs."""
from argparse import ArgumentParser
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import json,subprocess,sys
from flint import arb,ctx
D=Path(__file__).resolve().parent
p=ArgumentParser(description=__doc__)
for n in ('computation','tree','table32','table64','work'):p.add_argument(n,type=Path)
p.add_argument('--workers',type=int,default=2)
a=p.parse_args();C=a.computation.resolve();T=a.tree.resolve();A=a.table32.resolve();B=a.table64.resolve();W=a.work.resolve();W.mkdir(parents=True,exist_ok=False);assert a.workers>0

def run(script,*args):
 cmd=[sys.executable,'-E','-B',str(D/script),*map(str,args)];print('RUN',script,flush=True);subprocess.run(cmd,check=True)
def read(name):return json.loads((W/name).read_bytes())
L=F('.99999900555134042347742929191955');U=F('.999999005551341');G=F(120689264254125501,2**57)
run('all_materials.py',C,T,A,B,str(L),W/'materials.json');mat=read('materials.json');assert F(mat['scale'])==L and len(mat['results'])==10 and sorted((r['table'],r['order']) for r in mat['results'] if not r['positive'])==[(32,3),(64,3)]
run('triangle_threshold.py',C,T,A,B,W/'triangle.json');tri=read('triangle.json');assert F(tri['lower_s'])==L and F(tri['upper_s'])==U and len(tri['results'])==2
run('radial_binding.py',C,A,W/'radial.json');rad=read('radial.json');assert rad['zero_angle_and_counts_equal'] and rad['source_table_sha256']==mat['inputs']['3']
run('dense_pair_replay.py',C,T,B,W/'dense.json',str(a.workers));dense=read('dense.json');assert dense['all_passed'] and F(dense['lower_s'])==L and len(dense['leaves'])==24
run('cap_final.py',C,W/'cap.json');cap=read('cap.json');assert cap['complete'] and F(cap['scale'])==U and F(cap['projection_rate'])==F(57663128211,124250000000)
assert (W/'cap.json').read_bytes()==(D/'cap-final.json').read_bytes()
run('displacements_final.py',C,W/'points');pts=json.loads((W/'points/proof.json').read_bytes());assert pts['complete'] and pts['angular_rows']==2493 and pts['all_pairs']==435
run('interval_uniform.py',C,T,W/'interval.json',W/'points/proof.json');interval=read('interval.json');assert interval['complete'] and F(interval['lower'])==L and F(interval['upper'])==U and interval['native_leaves']==3273
ctx.prec=256;aa=arb(81827130)/2**28;cc=aa*(1+1/arb(3).sqrt());s0=3*(cc+(cc*cc+arb(16)*aa*aa/3).sqrt())/4;lower=arb(G.numerator)/G.denominator/s0**2;floor=F('.837451665604544468510297997471');assert lower>arb(floor.numerator)/floor.denominator
out={'complete':True,'scope':'conditional on frozen .83745 original-area trees and their geometric reduction','s_lower':str(L),'s_upper':str(U),'s0':str(s0),'goal':str(G),'bound':str(lower),'strict_decimal_floor':'0.837451665604544468510297997471','native_leaves':3273,'dense_leaves':24,'angular_rows':2493,'point_pairs':435,'hashes':{str(x.relative_to(W)):sha256(x.read_bytes()).hexdigest() for x in W.rglob('*.json')}}
(W/'result.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
