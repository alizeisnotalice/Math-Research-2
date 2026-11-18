"""Replay the deterministic P30 tree with an optional paired affine cap screen."""
from argparse import ArgumentParser
import ctypes
from fractions import Fraction as F
import json
from pathlib import Path
import subprocess
import sys
from time import monotonic

import numpy as np
from flint import arb

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
from certify_paired_affine_cap import KEY,GOAL,digest,fan_source,replace,compile_drivers

WRAPPER=r"""
extern "C" int paired_tree(const L* lower,const L* upper,const L* original,
    const L* anchor_lower,const L* anchor_upper,const std::int32_t* table,
    const std::int32_t* counts,L goal,int depth,int cap,const L* dirs,
    const L* disk,const L* endpoints,const L* pointers,const L* triples,
    int normals,int cells,L* out){
    Tester tester(table,counts,original,goal);
    Box lo,hi,whole_lo,whole_hi;
    std::copy(lower,lower+11,lo.begin());std::copy(upper,upper+11,hi.begin());
    std::copy(anchor_lower,anchor_lower+11,whole_lo.begin());
    std::copy(anchor_upper,anchor_upper+11,whole_hi.begin());
    for(int i=0;i<11;++i)
        if(whole_lo[i]>lo[i] || lo[i]>hi[i] || hi[i]>whole_hi[i])return 1;
    DynamicGuard run{tester,dirs,disk,endpoints,pointers,triples,whole_lo,whole_hi,
                     normals,cells,depth,cap,2};
    run.enable_pairing=true;
    run.init_collective();
    bool closed=run.walk(lo,hi,0);
    out[0]=closed;out[1]=run.calls;out[2]=run.old;out[3]=run.added;
    out[4]=run.depth;out[5]=run.reason;out[6]=run.tried;
    out[7]=run.paired_boxes;out[8]=run.paired_tries;out[9]=run.paired_cells;
    out[10]=run.paired_rows;out[11]=run.paired_tangents;
    if(closed && out[1]!=2*(out[2]+out[3])-1)std::abort();
    if(!closed){
        std::copy(run.unknown_lo.begin(),run.unknown_lo.end(),out+12);
        std::copy(run.unknown_hi.begin(),run.unknown_hi.end(),out+23);
    }
    return 0;
}
"""


def paired_driver(comp,work,cache):
    source=comp/'tomar_six_083730_guard.cpp'
    assert digest(source)=='42a4ee583d22940b3cd0b05f04285ba7e3d351745de707476211dce4efba5f81'
    source=source.read_text()
    source=replace(source,'#include "fan-exposed.hpp"','#include "paired-fan.hpp"')
    extra=(HERE/'tomar_paired_affine_tree.inc').read_text()
    source=replace(source,'    bool excluded(const Box& lo,const Box& hi){',
                   extra+'\n    bool excluded(const Box& lo,const Box& hi){')
    source=replace(source,'        bool needs_collective=false;',
                   '        bool needs_collective=false, used_pairing=false;')
    source=replace(source,'            if(!ok)return false;',r"""
            if(!ok){
                if(!enable_pairing)return false;
                if(!used_pairing){
                    ++paired_tries;
                    pair_lo=lo;pair_hi=hi;prepare_pairing();
                    if(min_plane(affine_area)>t.goal){++paired_boxes;return true;}
                }
                used_pairing=true;
                if(!paired_cell(i,L(e),support_upper,threshold))return false;
            }""")
    source=replace(source,'        if(needs_collective)++collective_boxes;',
              '        if(needs_collective)++collective_boxes;\n        if(used_pairing)++paired_boxes;')
    generated=work/'paired-tree.cpp'
    generated.write_text(source+WRAPPER)
    library=work/'paired-tree.so'
    subprocess.run(['clang++','-O3','-std=c++17','-ftrapv','-shared','-fPIC',
                    f'-fmodules-cache-path={cache}','-I',str(work),str(generated),
                    '-o',str(library)],check=True)
    driver=ctypes.CDLL(str(library)).paired_tree
    vec=lambda dtype:np.ctypeslib.ndpointer(dtype=dtype,ndim=1,flags='C_CONTIGUOUS')
    I,S,L,C=vec(np.int64),vec(np.int32),ctypes.c_int64,ctypes.c_int
    driver.argtypes=[I,I,I,I,I,S,S,L,C,C,I,I,I,I,I,C,C,I]
    driver.restype=C
    return driver,digest(generated)


def run(args):
    sys.path.insert(0,str(args.comp.resolve()))
    from certify_tomar_six_high import FINE,Q,SourceTree,build_kernel,load_table,ceil
    from replay_finite_point_083732 import PointProof
    from certify_nine import ar

    read=lambda name:json.loads((args.comp/name).read_bytes())
    data=read('finite-point-lift-083742.json')
    ledger=args.comp/'finite-point-lift-083742.json'
    assert digest(ledger)=='0e84a6fe33567f0784f08b88db7cae6dc79ef6807891bc175d510900139c8cec'
    for name,key in (('tomar-six-high-input.json','base_sha256'),
                     ('finite-point-leaf-input.json','point_sha256'),
                     ('finite-point-global-input.json','bank_sha256')):
        assert digest(args.comp/name)==data[key]
    base=read('tomar-six-high-input.json')
    tree=SourceTree(args.tree)
    table=load_table(args.table64,tree,base['tables']['64'])
    cloud=read('finite-point-leaf-input.json')
    points=set(map(tuple,cloud['points']))
    assert len(points)==30 and points=={(x,-y) for x,y in points}
    bank=read('finite-point-global-input.json')
    proof=PointProof(cloud,bank['cells'],100000000008865717121)
    pi=arb.pi()
    radius=ar(F(707,1000))
    assert (radius*radius<ar(F(1,2)) and
            (pi-(1/(2*radius)).acos())/4+(radius*radius-ar(F(1,4))).sqrt()/2>ar(F('.838')))

    args.work.mkdir(parents=True,exist_ok=True)
    fan,old,fan_sha=compile_drivers(args.comp,args.work,args.cache,base,data,build_kernel)
    paired,paired_sha=paired_driver(args.comp,args.work,args.cache)
    rootbox=tree.leaf(*KEY)
    support=proof.presolve(rootbox)
    normals,cells=len(support[1]),len(support[3])-1
    assert cells==384
    origin=np.array([int(F(x)*FINE) for x in tree.root['hi']],np.int64)
    A,B=[np.array([int(x*FINE) for x in row],np.int64) for row in zip(*rootbox,strict=True)]
    lo,hi=A.copy(),B.copy()
    guard_path=read('paired-affine-cap-low-depth44-input.json')['path']
    assert 0<=args.path_depth<=len(guard_path)
    path=guard_path[:args.path_depth]
    tracking=ceil(F('.83742')*2*Q*Q)

    def search(kernel,lo,hi,goal,depth,limit):
        out=np.zeros(34,np.int64)
        assert not kernel(lo,hi,origin,A,B,*table,goal,depth,limit,*support,normals,cells,out)
        return out

    for code in path:
        at=search(old,lo,hi,tracking,0,1)
        assert not at[0] and at[5]==2
        coefficients=np.zeros(25,np.int64)
        assert not fan(lo,hi,origin,*table,tracking,coefficients)
        i,side=divmod(code,2)
        assert i==coefficients[2] and lo[i]<hi[i] and (lo[i]+hi[i])%2==0
        middle=(lo[i]+hi[i])//2
        if side:lo[i]=middle
        else:hi[i]=middle

    cert=dict(
        leaf=list(KEY), path=path, box=[lo.tolist(),hi.tolist()], target=str(GOAL),
        depth_limit=args.depth, node_limit=args.limit, angular_cells=cells,
        tree_sha256=digest(ledger), paired_source_sha256=paired_sha, fan_sha256=fan_sha,
        driver_sha256=digest(Path(__file__))
    )
    meas={}
    goal=ceil(GOAL*2*Q*Q)
    for tag,kernel in (('original',old),('paired',paired)):
        if tag=='original' and args.paired_only:
            continue
        start=monotonic()
        print(f'{tag} begins at path depth {len(path)}',file=sys.stderr,flush=True)
        result=search(kernel,lo,hi,goal,args.depth,args.limit)
        meas[tag]=monotonic()-start
        cert[tag]=result[:7].tolist()
        print(f'{tag} ends after {meas[tag]:.3f}s: {cert[tag]}',file=sys.stderr,flush=True)
        if result[0]:
            assert result[1]==2*(result[2]+result[3])-1
        if tag=='paired':
            cert['paired_counters']=result[7:12].tolist()
            if not result[0]:
                cert['unresolved_box']=result[12:].tolist()
    if args.times:
        args.times.write_text(json.dumps(meas,separators=(',',':'))+'\n')
    return cert


if __name__=='__main__':
    assert __debug__
    parser=ArgumentParser(description=__doc__)
    for name in ('comp','tree','table64','work','cache'):
        parser.add_argument(name,type=Path)
    parser.add_argument('--path-depth',type=int,default=0)
    parser.add_argument('--depth',type=int,default=135)
    parser.add_argument('--limit',type=int,default=6_500_000)
    parser.add_argument('--times',type=Path)
    parser.add_argument('--paired-only',action='store_true')
    args=parser.parse_args()
    assert 1<=args.limit<=2**26 and 1<=args.depth<=135
    print(json.dumps(run(args),separators=(',',':')))
