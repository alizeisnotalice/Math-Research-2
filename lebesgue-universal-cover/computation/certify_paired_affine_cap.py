"""Verify a paired affine area/support cap on a true Tomar/P30 tree node."""
from argparse import ArgumentParser
import ctypes
from fractions import Fraction as F
from hashlib import sha256
from itertools import product
import json
from math import isqrt
from pathlib import Path
import subprocess
import sys

import numpy as np
from flint import arb

KEY=(1202,32259)
GOAL=F('.83743')
SLOPE_SCALE=10**10
BOUND_SCALE=10**12
DIM=11
FAN_WRAPPER=r"""
#include "paired-fan.hpp"

extern "C" int paired_read_fan(const L* lower,const L* upper,const L* original,
    const std::int32_t* table,const std::int32_t* counts,L goal,L* output){
    Tester tester(table,counts,original,goal);
    Box lo,hi;
    std::copy(lower,lower+11,lo.begin());std::copy(upper,upper+11,hi.begin());
    if(tester.check(lo,hi)<0)return 1;
    output[0]=tester.last;
    output[1]=tester.affine_center;
    output[2]=tester.split_axis;
    for(int i=0;i<11;++i){
        output[3+i]=tester.affine_anchor[i];
        output[14+i]=tester.affine_coeff[i];
    }
    return 0;
}
"""


BASE_WRAPPER=r"""
#include "tomar_six_083730_guard.cpp"

extern "C" int paired_run_original(const L* lower,const L* upper,const L* original,
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
    run.init_collective();
    bool closed=run.walk(lo,hi,0);
    out[0]=closed;out[1]=run.calls;out[2]=run.old;out[3]=run.added;
    out[4]=run.depth;out[5]=run.reason;out[6]=run.tried;
    return 0;
}
"""


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def replace(text, token, value):
    assert text.count(token)==1
    return text.replace(token,value)


def fan_source(text):
    text=replace(text,'W last=0; L goal;',
        'W last=0; L goal,affine_center=0,affine_coeff[11]{},affine_anchor[11]{};')
    text=replace(text,'last=0;pts.clear();',
        'last=0;affine_center=0; for(int i=0;i<11;++i)affine_coeff[i]=affine_anchor[i]=0; pts.clear();')
    text=replace(text,'ty=(lo[idx+1]+hi[idx+1])/(2*ratio);',
        '''ty=(lo[idx+1]+hi[idx+1])/(2*ratio);
        affine_anchor[idx]=tx*ratio*1024; affine_anchor[idx+1]=ty*ratio*1024;''')
    text=replace(text,'L anchor=L(row)*original_hi[k];',
        'L anchor=L(row)*original_hi[k]; affine_anchor[k]=anchor;')
    text=replace(text,'for(int k=1;k<5;++k){\n            error+=W(std::abs(g[k].y))*hx[k]+W(std::abs(g[k].x))*hy[k];','''
        W linear_losses=0;
        auto exact64=[](W v){
            if(v>std::numeric_limits<L>::max() || v<std::numeric_limits<L>::min())std::abort();
            return L(v);
        };
        for(int k=1;k<5;++k){
            int idx=(k-1)*3;
            affine_coeff[idx]=exact64(W(g[k].y)*Q);
            affine_coeff[idx+1]=exact64(-W(g[k].x)*Q);
            if(k>1) affine_coeff[idx-1]=exact64(grot[k]);
            W linear_error=W(std::abs(g[k].y))*hx[k]+W(std::abs(g[k].x))*hy[k];
            linear_losses+=linear_error;
            error+=linear_error;
    ''')
    text=replace(text,'twice+=scalar>=0?scalar/Q:-up(-scalar,Q);','''
            W linear_min=scalar>=0?scalar/Q:-up(-scalar,Q);
            linear_losses-=linear_min;
            twice+=linear_min;
    ''')
    return replace(text,'last=twice-error;',
        'last=twice-error;affine_center=exact64(last+linear_losses);')


def compile_drivers(comp, work, cache, base, receipt, kernel):
    work.mkdir(parents=True,exist_ok=True)
    kernel(work,base,cache)
    fan=work/'paired-fan.hpp'
    fan.write_text(fan_source((work/'fan-exposed.hpp').read_text()))
    wrapper=work/'paired-fan.cpp'
    wrapper.write_text(FAN_WRAPPER)
    original=comp/'tomar_six_083730_guard.cpp'
    assert digest(original)==receipt['weighted_sha256']
    base_wrapper=work/'paired-original.cpp'
    base_wrapper.write_text(BASE_WRAPPER)
    flags=['clang++','-O3','-std=c++17','-ftrapv','-shared','-fPIC',f'-fmodules-cache-path={cache}',
           '-I',str(work),'-I',str(comp.resolve())]
    for src,file in ((wrapper,'paired.so'),(base_wrapper,'original.so')):
        subprocess.run([*flags,str(src),'-o',str(work/file)],check=True)
    vector=lambda dtype:np.ctypeslib.ndpointer(dtype=dtype,ndim=1,flags='C_CONTIGUOUS')
    I,S,L,C=vector(np.int64),vector(np.int32),ctypes.c_int64,ctypes.c_int
    affine=ctypes.CDLL(str(work/'paired.so')).paired_read_fan
    affine.argtypes=[I,I,I,S,S,L,I]
    affine.restype=C
    original=ctypes.CDLL(str(work/'original.so')).paired_run_original
    original.argtypes=[I,I,I,I,I,S,S,L,C,C,I,I,I,I,I,C,C,I]
    original.restype=C
    return affine,original,digest(fan)


def min_affine(row):
    return row[0]+sum(min(x,0) for x in row[1:])


def max_affine(row):
    return row[0]+sum(max(x,0) for x in row[1:])


def affine_sum(rows):
    return tuple(map(sum,zip(*rows,strict=True)))


def scale(row,coefficient):
    return tuple(coefficient*x for x in row)


def sub(a,b):
    return tuple(x-y for x,y in zip(a,b,strict=True))


def value(row,point):
    return row[0]+sum(x*y for x,y in zip(row[1:],point,strict=True))


def largest_supports(values):
    keep=[]
    for i,a in enumerate(values):
        if all(i==j or (a==b and i<j) or (a!=b and min_affine(sub(b,a))<0)
               for j,b in enumerate(values)):
            keep.append(a)
    assert keep
    return keep


def unit_lengths(n2):
    num=n2.numerator*BOUND_SCALE**2
    down=isqrt(num//n2.denominator)
    return F(down,BOUND_SCALE),F(down+1,BOUND_SCALE)


def upper_depth(wall, eta):
    """Smallest upward rounded depth satisfying the one-wall sufficient test."""
    low,high=0,BOUND_SCALE//4
    def sufficient(d):
        depth=F(d,BOUND_SCALE)
        return depth**4>=4*eta**2*((wall+depth)**2-F(1,4))
    assert wall>=F(1,2) and eta>0 and sufficient(high)
    while high-low>1:
        mid=(low+high)//2
        if sufficient(mid):
            high=mid
        else:
            low=mid
    return F(high,BOUND_SCALE)


def fourth_root_upper(x):
    num=x.numerator*BOUND_SCALE**4
    down=isqrt(isqrt(num//x.denominator))
    return F(down+1,BOUND_SCALE)


def old_support_affines(direction, disk, endpoints, rootbox, lower, width, Q):
    parts=[(F(disk,Q),)+(F(0),)*DIM]
    for i in range(4):
        j=3*i
        z=[F(0)]*DIM
        z[j:j+2]=direction
        intercept=F(endpoints[2*i],Q)
        if i:
            j-=1
            left,right=rootbox[j]
            z[j]=F(endpoints[2*i+1]-endpoints[2*i],Q)/(right-left)
            intercept-=z[j]*left
        parts.append((intercept+sum(x*y for x,y in zip(z,lower,strict=True)),
                      *(x*y for x,y in zip(z,width,strict=True))))
    return largest_supports(parts)


def verify(args):
    sys.path.insert(0,str(args.comp.resolve()))
    from certify_tomar_six_high import FINE,Q,U,SourceTree,build_kernel,load_table,ceil
    from replay_finite_point_083732 import PointProof
    from certify_nine import ar
    record=json.loads(args.input.read_bytes())
    assert record['leaf']==list(KEY) and F(record['target'])==GOAL
    proof_goal=F(record['proof_goal'])
    assert GOAL<=proof_goal<F('.838')
    main=args.comp
    read=lambda name:json.loads((main/name).read_bytes())
    header_path=main/'finite-point-lift-083742.json'
    assert digest(header_path)=='0e84a6fe33567f0784f08b88db7cae6dc79ef6807891bc175d510900139c8cec'
    header=read(header_path.name)
    assert F(header['target'])==F('.83742')
    expected=next(r for r in header['rows'] if tuple(r[:2])==KEY)
    assert expected[2:]==[64,1360089,32642,647403,67,1316653,647403]
    for name,key in (('tomar-six-high-input.json','base_sha256'),
                     ('finite-point-leaf-input.json','point_sha256'),
                     ('finite-point-global-input.json','bank_sha256')):
        assert digest(main/name)==header[key]
    base=read('tomar-six-high-input.json')
    tree=SourceTree(args.tree)
    table=load_table(args.table64,tree,base['tables']['64'])
    cloud=read('finite-point-leaf-input.json')
    points=set(map(tuple,cloud['points']))
    assert len(points)==30 and points=={(x,-y) for x,y in points}
    bank=read('finite-point-global-input.json')
    point_proof=PointProof(cloud,bank['cells'],100000000008865717121)
    pi=arb.pi()
    radius=ar(F(707,1000))
    assert radius*radius<ar(F(1,2))
    assert ((pi-(1/(2*radius)).acos())/4+(radius*radius-ar(F(1,4))).sqrt()/2
            >ar(F('.838')))
    rootbox=tree.leaf(*KEY)
    support=point_proof.presolve(rootbox)
    scale_q=2*Q*Q
    goal=ceil(GOAL*scale_q)
    tracking_goal=ceil(F('.83742')*scale_q)
    main_upper=np.array([int(F(x)*FINE) for x in tree.root['hi']],np.int64)
    affine,baseline,fan_sha=compile_drivers(main,args.work,args.cache,base,header,build_kernel)

    anchor_lo,anchor_hi=[np.array([int(x*FINE) for x in row],np.int64)
                          for row in zip(*rootbox,strict=True)]

    def run_old(lo,hi,target,depth):
        out=np.zeros(30,np.int64)
        assert not baseline(lo,hi,main_upper,anchor_lo,anchor_hi,*table,target,depth,
                            2**(depth+1)-1,*support,len(support[1]),len(support[3])-1,out)
        return out

    def run_fan(lo,hi,target):
        out=np.zeros(25,np.int64)
        assert not affine(lo,hi,main_upper,*table,target,out)
        return list(map(int,out))

    lo,hi=anchor_lo.copy(),anchor_hi.copy()
    for code in record['path']:
        old=run_old(lo,hi,tracking_goal,0)
        assert old[0]==0 and old[5]==2
        axis,side=divmod(code,2)
        assert run_fan(lo,hi,tracking_goal)[2]==axis
        assert (lo[axis]+hi[axis])%2==0 and lo[axis]<hi[axis]
        mid=(lo[axis]+hi[axis])//2
        if side:
            lo[axis]=mid
        else:
            hi[axis]=mid
    assert [lo.tolist(),hi.tolist()]==record['box']
    parent=run_old(lo,hi,goal,0)
    old_expected=record['original_subtree']
    assert old_expected[0]==1 and old_expected[1]>1 and 1<=old_expected[4]<=20 and old_expected[5]==0
    children=run_old(lo,hi,goal,old_expected[4])
    assert parent[0]==0 and parent[1]==1 and parent[5]==2
    assert children[:7].tolist()==old_expected, children[:7].tolist()
    assert int(children[1])==2*int(children[2]+children[3])-1
    low_five=None
    if record.get('low_five_pose', False):
        from certify_five_harmonic import placement_area
        pose_file=main/'mixed-support-blind-low-pose.json'
        pose=json.loads(pose_file.read_bytes())
        assert pose['coordinate_denominator']==FINE
        z=pose['pose']
        assert all(l<=x<=h for l,x,h in zip(lo,z,hi,strict=True))
        area_upper=F(placement_area([F(x,FINE) for x in z],pi,k=1024)['upper'])
        assert area_upper<F('.837345331409')<GOAL
        low_five=dict(five_area_upper=str(F(ceil(area_upper*BOUND_SCALE),BOUND_SCALE)),
                      pose_numerator=z,pose_sha256=digest(pose_file))
    fan=run_fan(lo,hi,goal)
    lower=[F(int(v),FINE) for v in lo]
    width=[F(int(v),FINE) for v in hi-lo]
    coeff=[F(v,scale_q) for v in fan[14:25]]
    anchor=[F(v,1024*FINE) for v in fan[3:14]]
    area=(F(fan[1],scale_q)+sum(z*(x-y) for z,x,y in zip(coeff,lower,anchor,strict=True)),
          *(z*x for z,x in zip(coeff,width,strict=True)))
    floor=min_affine(area)
    eta=proof_goal-min(floor,F(fan[0],scale_q))
    raw_directions=np.asarray(support[0]).reshape(-1,2).tolist()
    anchors=np.asarray(support[2]).reshape(8,-1).T.tolist()
    n2=[F(x*x+y*y,U*U) for x,y in raw_directions]
    lengths=[unit_lengths(z) for z in n2]
    affines=[old_support_affines((F(n[0],U),F(n[1],U)),int(d),h,rootbox,lower,width,Q)
             for n,d,h in zip(raw_directions,support[1],anchors,strict=True)]
    maxes=[max(map(max_affine,x)) for x in affines]
    walls=[b/l[0] for b,l in zip(maxes,lengths,strict=True)]
    depths=[upper_depth(b,eta) for b in walls]
    scale_i=[ell[1]*fourth_root_upper(4*((b+d)**2-F(1,4)))
             for ell,b,d in zip(lengths,walls,depths,strict=True)]
    pointers=list(map(int,support[3]))
    rows=np.asarray(support[4]).reshape(-1,8).tolist()
    assert len(record['rows'])==len(pointers)-1==384
    all_bounds=[]
    separate=[]
    comparison=None
    for cell,(pick,lm) in enumerate(record['rows']):
        assert type(pick) is type(lm) is int and pointers[cell]<=pick<pointers[cell+1] and lm>=0
        chosen=None
        independent=[]
        # A simple independent cap bound uses the same affine fan,
        # downward support endpoints, and upper quadratic costs.
        for rid in range(pointers[cell],pointers[cell+1]):
            r=rows[rid]
            inds=r[:3]
            weights=[F(x,r[6]) for x in r[3:6]]
            ns=[raw_directions[j] for j in inds]
            assert sum(weights)==1 and min(weights)>0
            assert all(sum(w*n[i] for w,n in zip(weights,ns,strict=True))==0 for i in (0,1))
            factors=[w*scale_i[j] for j,w in zip(inds,weights,strict=True)]
            acute=[(i,j) for i,j in ((0,1),(0,2),(1,2))
                   if ns[i][0]*ns[j][0]+ns[i][1]*ns[j][1]>0]
            assert len(acute)<=1
            C=sum(f*f for f in factors)+sum(2*factors[i]*factors[j] for i,j in acute)
            H=F(r[7],Q)
            m=H-sum(w*maxes[j] for j,w in zip(inds,weights,strict=True))
            independent.append(floor+max(m,0)**2/C)
            if rid==pick:
                chosen=(inds,weights,C,H)
        inds,weights,C,H=chosen
        slope=F(lm,SLOPE_SCALE)
        base_row=(area[0]+2*slope*H-slope*slope*C,*area[1:])
        terms=[tuple(scale(row,-2*slope*w) for row in affines[j])
               for j,w in zip(inds,weights,strict=True)]
        combined=min(min_affine(affine_sum((base_row,*items))) for items in product(*terms))
        assert combined>proof_goal, (cell, float(combined), float(proof_goal))
        all_bounds.append(combined)
        separate.append(max(independent))
        if separate[-1]<min(separate[:-1],default=F(1)) and separate[-1]<GOAL:
            a=tuple(int(x<0) for x in area[1:])
            options=[tuple(scale(row,w) for row in affines[j])
                     for j,w in zip(inds,weights,strict=True)]
            top=max((affine_sum(x) for x in product(*options)),key=max_affine)
            b=tuple(int(x>0) for x in top[1:])
            assert a!=b
            def describe(point):
                h=sum(w*max(value(row,point) for row in affines[j])
                      for j,w in zip(inds,weights,strict=True))
                cap=2*slope*(H-h)-slope*slope*C
                return {
                    'position_numerator': [int(v) for v in lo+np.asarray(point)*(hi-lo)],
                    'five_affine': str(value(area,point)),
                    'cap_affine': str(cap),
                    'sum':str(value(area,point)+cap),
                }
            comparison=dict(cell=cell,row=pick,a=describe(a),b=describe(b),
                            independent=str(max(independent)),coupled_lower=str(combined))
    worst=min(all_bounds)
    assert min(separate)<GOAL and comparison is not None
    return dict(
        leaf=list(KEY),path=record['path'],box=record['box'],target=str(GOAL),proof_goal=str(proof_goal),
        reference_tree_sha256=digest(header_path),
        original_parent=parent[:7].tolist(),original_full_subtree=children[:7].tolist(),
        original_split_axis=fan[2],paired_nodes=1,p30_angle_cells=len(all_bounds),
        coupled_conditional_floor_at_1e12=str(F((worst*BOUND_SCALE).__floor__(),BOUND_SCALE)),
        uncoupled_min_at_1e12=str(F((min(separate)*BOUND_SCALE).__ceil__(),BOUND_SCALE)),
        weakest_cell=all_bounds.index(worst),
        witness_pair=comparison,fan_floor=str(floor),source_sha256=digest(Path(__file__)),
        paired_fan_sha256=fan_sha,input_sha256=digest(args.input),certified_low_five=low_five
    )


if __name__=='__main__':
    assert __debug__
    parser=ArgumentParser(description=__doc__)
    for name in ('comp','tree','table64','work','cache','input'):
        parser.add_argument(name,type=Path)
    print(json.dumps(verify(parser.parse_args()),separators=(',',':')))
