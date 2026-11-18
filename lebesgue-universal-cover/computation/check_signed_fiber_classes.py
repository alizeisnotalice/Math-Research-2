"""Screen one fixed original leaf by a dominant translation with exact signed area floors."""
from argparse import ArgumentParser
import ctypes
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
import subprocess
import sys

import numpy as np

TARGET=F('.83743')
PART=1202
POS=32259
AREA=2**57
FINE=2**44
CPP=r"""
void screen_box(Tester& t,Box& lo,Box& hi,L* out,L* boxes,int depth,int maxdepth) {
    ++out[0];
    const int ok=t.check(lo,hi);
    if(ok<0)std::abort();
    if(ok>0) {++out[1];if(t.old_observed<t.goal)++out[36];return;}
    const int axis=t.split_axis;
    const L a=lo[axis],b=hi[axis];
    if(depth<maxdepth && b-a>=2) {
        const L mid=a+(b-a)/2;
        hi[axis]=mid;screen_box(t,lo,hi,out,boxes,depth+1,maxdepth);hi[axis]=b;
        lo[axis]=mid;screen_box(t,lo,hi,out,boxes,depth+1,maxdepth);lo[axis]=a;
        return;
    }
    const L at=out[2];
    std::copy(lo.begin(),lo.end(),boxes+22*at);std::copy(hi.begin(),hi.end(),boxes+22*at+11);
    if(out[2]++==0) {
        std::copy(lo.begin(),lo.end(),out+3);std::copy(hi.begin(),hi.end(),out+14);
    }else{
        for(int j=0;j<11;++j){
            out[3+j]=std::min(out[3+j],lo[j]);
            out[14+j]=std::max(out[14+j],hi[j]);
        }
    }
}
extern "C" int classify(int n,const L* lower,const L* upper,const L* source_hi,
                  const std::int32_t* table,const std::int32_t* size,L goal,int depth,L* out,L* boxes){
    Tester tester(table,size,source_hi,goal);
    Box lo,hi;
    for(int i=0;i<n;++i){
        std::copy(lower+11*i,lower+11*i+11,lo.begin());
        std::copy(upper+11*i,upper+11*i+11,hi.begin());
        L* ptr=out+38*i;
        tester.check(lo,hi);
        ptr[25]=tester.observed;ptr[26]=tester.old_observed;ptr[27]=tester.split_axis;
        for(int j=0;j<8;++j)ptr[28+j]=tester.gradient[j];
        screen_box(tester,lo,hi,ptr,boxes+(L(22)*i<<depth),0,depth);
    }
    return 0;
}
"""


def hash_bytes(x):
    return sha256(x).hexdigest()


def compile_scan(comp,work,cache):
    import replay_tomar_083731_scored as reference
    from replay_tomar_signed_translation_screen import compile_screen

    work.mkdir(parents=True,exist_ok=True)
    def build(signed,work,cache):
        source=(signed/'five_grouped_localization.cpp').read_text().split('\nstruct Runner',1)[0]
        head='L goal;int total=0;int split_axis=0;'
        repl='L goal,observed=0,old_observed=0,gradient[8]{};int total=0;int split_axis=0;'
        assert source.count(head)==1
        source=source.replace(head,repl)
        line='int check(const Box& lo,const Box& hi){'
        assert source.count(line)==1
        source=source.replace(line,line+'\n        observed=old_observed=0;')
        token='// The signed area of this fixed material chain is separately affine'
        assert source.count(token)==1
        source=source.replace(token,'old_observed=small(twice-error);\n        '+token)
        end='return twice-error>=goal ? 1 : 0;'
        assert source.count(end)==1
        source=source.replace(end,'''observed=small(twice-error);
        for(int j=1;j<=4;++j){
           gradient[2*j-2]=g[j].y;
           gradient[2*j-1]=-g[j].x;
        }
        '''+end)
        raw=(source+'\n'+CPP).encode()
        code=work/'signed-fiber.cpp'; code.write_bytes(raw)
        library=work/'signed-fiber.dylib'
        subprocess.run(['clang++','-O3','-std=c++17','-ftrapv','-shared','-fPIC',
                        f'-fmodules-cache-path={cache}',str(code),'-o',str(library)],check=True)
        handle=ctypes.CDLL(str(library)).classify
        v=lambda d:np.ctypeslib.ndpointer(dtype=d,ndim=1,flags='C_CONTIGUOUS')
        A=v(np.int64);B=v(np.int32)
        handle.argtypes=[ctypes.c_int,A,A,A,B,B,ctypes.c_int64,ctypes.c_int,A,A]
        handle.restype=ctypes.c_int
        return handle,hash_bytes(raw)
    return compile_screen(reference,build,comp,work,cache)


def scan(a):
    sys.path.insert(0,str(a.comp))
    from certify_tomar_six_high import SourceTree
    from verify_tomar_fan_table import verified_table
    tree=SourceTree(a.tree)
    sample,shape=verified_table(a.table,a.tree,tree.root)
    material=np.ascontiguousarray(sample,dtype=np.int32).ravel()
    sizes=np.asarray(shape,dtype=np.int32)
    box=tree.leaf(PART,POS)
    lo,hi=(np.array([int(z*FINE) for z in x],np.int64) for x in zip(*box,strict=True))
    root_hi=np.array([int(F(x)*FINE) for x in tree.root['hi']],np.int64)
    f,driver=compile_scan(a.comp,a.work,a.cache)
    goal=int(TARGET*AREA)+1
    coarse=np.zeros(38,np.int64)
    assert not f(1,lo,hi,root_hi,material,sizes,goal,0,coarse,np.zeros(22,np.int64))
    translations=(0,1,3,4,6,7,9,10)
    influence=[int(abs(v))*int(hi[j]-lo[j]) for j,v in zip(translations,coarse[28:36],strict=True)]
    axis=translations[influence.index(max(influence))]
    assert int(hi[axis]-lo[axis])%a.bins==0
    boundaries=np.linspace(lo[axis],hi[axis],a.bins+1).astype(np.int64)
    step=(hi[axis]-lo[axis])//a.bins
    assert np.array_equal(boundaries, lo[axis]+step*np.arange(a.bins+1,dtype=np.int64))
    Ls=np.tile(lo,(a.bins,1));Hs=np.tile(hi,(a.bins,1))
    Ls[:,axis]=boundaries[:-1];Hs[:,axis]=boundaries[1:]
    out=np.zeros((a.bins,38),np.int64)
    shapes=np.zeros((a.bins,2**a.depth,22),np.int64)
    assert not f(a.bins,Ls.ravel(),Hs.ravel(),root_hi,material,sizes,goal,a.depth,out.ravel(),shapes.ravel())
    assert np.all(out[:,0]==2*(out[:,1]+out[:,2])-1)
    print(json.dumps(dict(key=[PART,POS],target=str(TARGET),grid=FINE,area_units=AREA,
                    root=[lo.tolist(),hi.tolist()],root_scores=coarse[25:28].tolist(),
                    signed_gradient_weight=influence,classification_axis=axis,
                    boundaries=boundaries.tolist(),depth=a.depth,slices=out[:,:28].tolist(),
                    new_signed_closures=out[:,36].tolist(),
                    live_boxes=[shapes[i,:int(o[2])].tolist() for i,o in enumerate(out)],
                    kernel_sha256=driver,source_sha256=hash_bytes(Path(__file__).read_bytes())),separators=(',',':')))


if __name__=='__main__':
    p=ArgumentParser(description=__doc__)
    for name in ('comp','tree','table','work','cache'):
        p.add_argument(name,type=Path)
    p.add_argument('--depth',type=int,default=11)
    p.add_argument('--bins',type=int,default=4)
    scan(p.parse_args())
