"""Certify complete dyadic parts of the original box (180,24225) at 0.83744."""
import argparse
import ctypes
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
import shutil
import subprocess
import sys

from flint import ctx
import numpy as np

TARGET=F('.83744')
GOAL=120687823102244742
PART,POS=180,24225


def main():
    if not __debug__:
        raise SystemExit('Run with assertions enabled.')
    ap=argparse.ArgumentParser(description=__doc__)
    for name in ('computation','tree','table64','work','cache','selection'):
        ap.add_argument(name,type=Path)
    args=ap.parse_args()
    sys.path.insert(0,str(Path(__file__).resolve().parent))
    comp=args.computation.resolve()
    sys.path.insert(0,str(comp))
    from certify_tomar_six_high import FINE,Q,SourceTree,build_kernel,ceil,load_table
    from replay_finite_point_083732 import PointProof
    ctx.prec=160
    assert ceil(2*Q*Q*TARGET)==GOAL and F(GOAL,2*Q*Q)>TARGET
    previous=(comp/'finite-point-lift-083742.json').read_bytes()
    assert sha256(previous).hexdigest()=='0e84a6fe33567f0784f08b88db7cae6dc79ef6807891bc175d510900139c8cec'
    old=json.loads(previous)
    def data(name,field):
        b=(comp/name).read_bytes();assert sha256(b).hexdigest()==old[field];return json.loads(b)
    base=data('tomar-six-high-input.json','base_sha256')
    bank=data('finite-point-global-input.json','bank_sha256')
    cloud=data('finite-point-leaf-input.json','point_sha256')
    proof=PointProof(cloud,bank['cells'],100000000008865717121)
    tree=SourceTree(args.tree)
    mat=load_table(args.table64,tree,base['tables']['64'])
    def scaled(xx):
        q=[F(v)*FINE for v in xx];assert all(n.denominator==1 for n in q)
        return np.array([int(n)for n in q],dtype=np.int64)
    lo,hi=(scaled(k) for k in zip(*tree.leaf(PART,POS),strict=True))
    original=scaled(tree.root['hi'])
    work=args.work.resolve();cache=args.cache.resolve()
    build_kernel(work,base,cache)
    src=comp/'tomar_six_083730_guard.cpp'
    assert sha256(src.read_bytes()).hexdigest()==old['weighted_sha256']
    so=work/'prefix-guard.so'
    subprocess.run([shutil.which('clang++'),'-O3','-std=c++17','-ftrapv','-shared','-fPIC',
                    f'-fmodules-cache-path={cache}','-I',str(work),src,'-o',so],check=True)
    run=ctypes.CDLL(str(so)).lift_guard_run
    v=lambda d:np.ctypeslib.ndpointer(dtype=d,ndim=1,flags='C_CONTIGUOUS')
    i64,i32=v(np.int64),v(np.int32)
    run.argtypes=[i64,i64,i64,i32,i32,ctypes.c_int64,ctypes.c_int,ctypes.c_int,
                  i64,i64,i64,i64,i64,ctypes.c_int,ctypes.c_int,ctypes.c_int,i64]
    run.restype=ctypes.c_int
    # The original arithmetical guard defines the inherited partition.
    # Only the terminal full-box test may use the separately checked mesh.
    mesh_lib=None
    def mesh_guard():
        nonlocal mesh_lib
        if mesh_lib is None:
            from tomar_p30_mesh_guard import write_kernel
            meshsrc=write_kernel(work,src)
            lib=work/'mesh-guard.so'
            subprocess.run([shutil.which('clang++'),'-O3','-std=c++17','-ftrapv','-shared','-fPIC',
                            f'-fmodules-cache-path={cache}','-I',str(work),meshsrc,'-o',lib],check=True)
            mesh_lib=ctypes.CDLL(str(lib)).lift_guard_run
            mesh_lib.argtypes=run.argtypes
            mesh_lib.restype=run.restype
        return mesh_lib
    def test(l,h,cap,depth,use_mesh=False,dense=False):
        materials=mat
        box=[(F(int(x),FINE),F(int(y),FINE)) for x,y in zip(l,h,strict=True)]
        if use_mesh:
            from tomar_p30_mesh_guard import piecewise_support,dense_material
            if dense:materials=dense_material(mat,l,h,original)
            from certify_tomar_six_high import U
            support=piecewise_support(proof,box,Q,U,ceil)
            guard=mesh_guard()
        else:
            support=proof.presolve(box)
            guard=run
        st=np.zeros(30,dtype=np.int64)
        assert guard(l,h,original,*materials,GOAL,depth,cap,*support,len(support[1]),len(support[3])-1,2,st)==0
        return st
    raw=args.selection.read_bytes();sel=json.loads(raw);paths=sel['paths']
    assert sel['part']==PART and sel['pos']==POS and sel['target']=='.83744'
    assert paths and len(paths)==len(set(paths)) and all(set(s)<={'0','1'} and len(s)<32 for s in paths)
    variant=sel.get('terminal_guard','')
    use_mesh=variant in ('mesh64','mesh64-dense')
    rows=[]
    for bits in paths:
        l,h=lo.copy(),hi.copy()
        steps=[]
        for bit in bits:
            st=test(l,h,1,155)
            if st[0]:
                break
            assert st[5]==1 and st[1]==1
            left_lo,left_hi=st[7:18],st[18:29]
            axes=np.flatnonzero(left_hi!=h)
            assert np.array_equal(left_lo,l) and len(axes)==1
            k=int(axes[0]);m=int(left_hi[k])
            assert int(l[k])<m<int(h[k]) and 2*m==int(l[k])+int(h[k])
            if len(steps)>=10:
                # P30 interpolation uses only these three rotation intervals.
                # Halving the widest one makes that bound sharper on both closed children.
                k=max((2,5,8),key=lambda j:int(h[j])-int(l[j]))
                total=int(l[k])+int(h[k])
                assert total%2==0
                m=total//2
                assert int(l[k])<m<int(h[k])
            steps.append([k,m])
            if bit=='0':h[k]=m
            else:l[k]=m
        else:
            st=test(l,h,268435455,240 if use_mesh else 155,use_mesh,variant=='mesh64-dense')
        if st[0]:
            assert st[1]==2*(st[2]+st[3])-1 and st[5]==0
        row=dict(bits=bits,steps=steps,lo=l.tolist(),hi=h.tolist(),stat=st.tolist())
        if use_mesh:row['terminal_guard']=variant
        (work/f'prefix-{bits}.json').write_text(json.dumps(row,separators=(',',':'))+'\n')
        rows.append(row)
    result=dict(complete=True,all_closed=all(r['stat'][0]==1 for r in rows),part=PART,pos=POS,
                target=str(TARGET),threshold_integer=GOAL,original_lo=lo.tolist(),original_hi=hi.tolist(),
                selection_sha256=sha256(raw).hexdigest(),rows=rows)
    (work/'result.json').write_text(json.dumps(result,separators=(',',':'))+'\n')


if __name__=='__main__':
    main()
