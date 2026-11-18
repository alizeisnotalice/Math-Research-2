"""Independently build and verify integer five-body witnesses for Metal-proposed material moves."""
import argparse
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
import subprocess
import sys
import numpy as np
from flint import ctx

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
import tomar_dynamic_five as dynamic

def h(path):return sha256(path.read_bytes()).hexdigest()

def main():
    if not __debug__:raise SystemExit("Assertions required.")
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('upstream','cert','proposals','work','cache'):p.add_argument(key,type=Path)
    a=p.parse_args();a.work.mkdir(parents=True,exist_ok=True)
    assert h(a.cert/'part0000.json')=='e28aff763d77762b968bcc37a76c23dde326e6fbf6923589a10cbd723229fee5'
    assert h(a.proposals)=='03dde654589a2f64974b193d066044fa9bd5409d5c71131a0fd882da04530ef8'
    revision='10af98606c98ca9ffd0d9f402cae1e2efe74b792'
    for name in ['gen8.py','gen11.py']:
        assert (a.upstream/'independent'/name).read_bytes()==subprocess.check_output(['git','-C',a.upstream,'show',revision+':independent/'+name])
    sys.path.insert(0,str(a.upstream.resolve()/'independent'));import gen11
    lib=a.work/'five-affine.so'
    subprocess.run(['clang++','-O3','-std=c++17','-ftrapv',f'-fmodules-cache-path={a.cache}','-shared','-fPIC',
                    HERE/'tomar_dynamic_affine.cpp','-o',lib],check=True)
    target=Fraction('0.83742');ctx.prec=144;dynamic.configure(a.upstream,a.cert,lib,target)
    assert dynamic.gen8.placements is gen11.placements
    from tomar_leaf_box import leaf_box
    raw=(a.cert/'SHA256SUMS').read_bytes()
    assert sha256(raw).hexdigest()=='56cbe5dfdadc6b5aa07cf6ec9acf4f637a3c4fa78db54df0ae2d78916d292410'
    manifests={n.removeprefix('./'):v for v,n in (line.split() for line in raw.decode().splitlines())}
    with np.load(a.proposals,allow_pickle=False) as z: proposal=z['proposals']
    assert sorted(set(proposal[:,2]))==list(range(23,31))
    active=None;rows=[];goal=dynamic.goal
    for part,idx,winner in proposal.tolist():
        if part!=active:
            stem=f'part{part:04}'
            for name in (stem+'.json',stem+'.npy'):
                assert h(a.cert/name)==manifests[name]
            header=json.loads((a.cert/(stem+'.json')).read_bytes())
            codes=np.load(a.cert/(stem+'.npy'),allow_pickle=False)
            assert part==header['part']
            prefix=np.r_[1,1+np.cumsum(np.where(codes<11,1,-1))]
            assert prefix[-1]==0 and prefix[:-1].min()==1
            active=part
        assert 0<=idx<len(codes) and codes[idx]>=253
        c,w=leaf_box(prefix,codes,part,idx,header)
        ends=[(Fraction(x)+s*Fraction(y))*2**44 for s in (-1,1) for x,y in zip(c,w)]
        assert all(e.denominator==1 for e in ends)
        lo,hi=np.asarray(ends,np.int64).reshape(2,11)
        pts,cnt=dynamic.table(lo,hi,True)
        def exact(move):
            out=np.zeros(12,np.int64)
            count=dynamic.dll.forms_for_box(lo,hi,dynamic.hi,pts,cnt,dynamic.moves[move],1,out)
            assert count in (0,1)
            return 0 if not count else int(out[0])-sum(abs(int(v))for v in out[1:])
        old=exact(0);new=exact(winner)
        rows.append([part,idx,winner,old,new])
    rows=np.array(rows,np.int64)
    out=a.work/'ranked-integer-floors.npz';np.savez_compressed(out,rows=rows)
    r=dict(threshold=str(target),strict_integer_floor=goal,candidates=len(proposal),center_only_pass=int((rows[:,3]>=goal).sum()),
      selected_move_pass=int((rows[:,4]>=goal).sum()),additional_vs_center=int(((rows[:,4]>=goal)&(rows[:,3]<goal)).sum()),
      best_of_center_and_selected=int((np.maximum(rows[:,3],rows[:,4])>=goal).sum()),
      selected_move_higher=int((rows[:,4]>rows[:,3]).sum()),upstream=revision,sha256=h(out),
      proposals_sha256=h(a.proposals),source_sha256=h(Path(__file__)))
    (a.work/'result.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r))
if __name__=='__main__':main()
