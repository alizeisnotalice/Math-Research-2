"""Recreate and compare old and lower-endpoint dense anchors on all 24 original leaves."""
from pathlib import Path
from fractions import Fraction as F
from concurrent.futures import ProcessPoolExecutor,as_completed
from hashlib import sha256
import sys,json
import numpy as np
D=Path(__file__).resolve().parent;C=Path(sys.argv[1]).resolve();TREE=Path(sys.argv[2]).resolve();TABLE=Path(sys.argv[3]).resolve();OUT=Path(sys.argv[4]).resolve();workers=int(sys.argv[5]) if len(sys.argv)>5 else 2
sys.path[:0]=[str(D),str(C),str(C/'six-input-083745')]
from certify_tomar_six_high import FINE,SourceTree,load_table
from replay_dense_083745 import inputs
from tomar_p30_scaled_material import scaled_material as old
from scaled_material_lower import scaled_material as new
oldpath=C/'tomar_p30_scaled_material.py';oldsrc=oldpath.read_text();assert sha256(oldpath.read_bytes()).hexdigest()=='da3fdc4472aaf0b09e020b6a3a7b467f227fe695c84ad0c4f502244e0f20cc2b';assert oldsrc.count('scale=arb(999999006)/1000000000')==1
assert (D/'scaled_material_lower.py').read_text()==oldsrc.replace('scale=arb(999999006)/1000000000','scale=arb(19999980111026808469548585838391)/20000000000000000000000000000000')
def job(item):
 key,raw,expected=item;data=json.loads(raw);tree=SourceTree(TREE);base=json.loads((C/'tomar-six-high-input.json').read_text());mat=load_table(TABLE,tree,base['tables']['64'])
 def scaled(xx):
  q=[F(v)*FINE for v in xx];assert all(n.denominator==1 for n in q);return np.array([int(n) for n in q],dtype=np.int64)
 lo,hi=(scaled(k) for k in zip(*tree.leaf(data['part'],data['pos']),strict=True));orig=scaled(tree.root['hi']);arr0,count0=old(mat,lo,hi,orig);arr1,count1=new(mat,lo,hi,orig);assert np.array_equal(count0,count1) and np.array_equal(arr0,arr1)
 return {'leaf':key,'selection_sha256':sha256(raw).hexdigest(),'dense_sha256':sha256(arr0.tobytes()).hexdigest(),'counts':count0.tolist(),'integer_points':int(len(arr0)//2),'original_receipt_closed':bool(expected['all_closed'])}
if __name__=='__main__':
 jobs=inputs(C/'six-input-083745/volumes');res=[]
 with ProcessPoolExecutor(max_workers=workers) as pool:
  for f in as_completed([pool.submit(job,x) for x in jobs]):
   x=f.result();res.append(x);print(x['leaf'],flush=True)
 assert len(res)==len({x['leaf'] for x in res})==24
 OUT.write_text(json.dumps({'lower_s':str(F('.99999900555134042347742929191955')),'all_passed':True,'leaves':sorted(res,key=lambda x:x['leaf'])},indent=2)+'\n')
