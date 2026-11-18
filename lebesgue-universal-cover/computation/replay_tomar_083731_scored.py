"""Localize the next five-body threshold and preserve exact integer leaf lower bounds."""
from argparse import ArgumentParser
from concurrent.futures import ThreadPoolExecutor
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import ctypes
import json
import os
import re
import subprocess
import sys
from time import monotonic

import numpy as np

from audit_tomar_outer_083730 import AREA, FINE, sha
from select_tomar_083731_scores import radius_rows

HERE=Path(__file__).resolve().parent
SOURCE_SHA='5db2c5462580c07b4d1e50092a41ca95df1ae455c560e7f0ae65d458f447691e'
TABLE_SHA='fe5b816060d713642dc44e3f8981f4ae5d1b2b1c1e271b7dff75f031dd75d87c'
STAGES=((2,127),(1,1023),(1,8191),(1,65535))


def kernel(main,work,cache):
    raw=(main/'five_grouped_localization.cpp').read_bytes()
    assert sha(raw)==SOURCE_SHA
    prefix=raw.decode().split('\nstruct Runner',1)[0]
    edits=[
        ('bool clip_domain(Box& lo,Box& hi){','bool clip_domain(Box& lo,Box& hi,const L* r){'),
        ('    constexpr L r[4]={1202425916136,1169880371954,1056234850107,1002578682672};\n',''),
        ('L goal;int total=0;int split_axis=0;', 'L goal,current_floor=0;int total=0;int split_axis=0;'),
        ('return twice-error>=goal ? 1 : 0;',
         'if(twice-error<goal)return 0;\n        current_floor=small(twice-error);return 1;'),
    ]
    for a,b in edits:
        assert prefix.count(a)==1
        prefix=prefix.replace(a,b)
    tail=(HERE/'tomar_083731_scored_tail.cpp').read_bytes()
    full=prefix.encode()+b'\n'+tail
    source=work/'outer-scored.cpp'
    lib=work/'outer-scored.dylib'
    source.write_bytes(full)
    compiler=os.environ.get('CXX','clang++')
    flags=['-O3','-std=c++17','-ftrapv','-shared','-fPIC',f'-fmodules-cache-path={cache}']
    subprocess.run([compiler,*flags,str(source),'-o',str(lib)],check=True)
    fn=ctypes.CDLL(str(lib)).adaptive_scored_trim
    v=lambda dtype:np.ctypeslib.ndpointer(dtype=dtype,ndim=1,flags='C_CONTIGUOUS')
    integer=ctypes.c_int
    l64=v(np.int64);l32=v(np.int32);u8=v(np.uint8)
    fn.argtypes=[u8,integer,l64,l64,l64,u8,u8,l64,l32,l32,ctypes.c_int64,integer,integer,
                 l64,l64,l64,l64]
    return fn,sha(full),sha(lib.read_bytes())


class Data:
    def __init__(self,main,tree,table_path,seed_file):
        from verify_tomar_fan_table import verified_table
        manifest=(tree/'SHA256SUMS').read_bytes()
        seed=json.loads(seed_file.read_bytes())
        self.seed=seed
        assert sha(manifest)==seed['tree_manifest_sha256']
        self.manifest={n.removeprefix('./'):h for h,n in (x.split() for x in manifest.decode().splitlines())}
        self.tree=tree
        self.seed_rows={p:r for p,*r in seed['rows']}
        self.header=json.loads(self.read('part0000.json'))
        self.lo,self.hi=(np.array([int(F(v)*FINE) for v in self.header[k]],np.int64) for k in ('lo','hi'))
        assert self.header['split0']==[2]*11
        assert np.all((self.lo+self.hi)%2==0)
        self.mid=(self.lo+self.hi)//2
        assert sha(table_path.read_bytes())==TABLE_SHA
        tab,counts=verified_table(table_path,tree,self.header)
        self.tables={}
        offsets=np.r_[0,np.cumsum(counts)]
        for div,_ in STAGES:
            ids=np.concatenate([np.arange(a,b,div) for a,b in zip(offsets[:-1],offsets[1:],strict=True)])
            self.tables[div]=(np.ascontiguousarray(tab[:,ids,:]).ravel(),
                              ((counts+div-1)//div).astype(np.int32))
        raw=(main/'five-body-radial-affine.json').read_bytes()
        assert sha(raw)==seed['provenance']['radial_certificate']
        rows=radius_rows(json.loads(raw),F(seed['radial_target']))
        assert rows==seed['radii']
        self.radii=np.array([r['integer_radius'] for r in rows],np.int64)
        self.radial_floors=np.array([r['integer_area_lower'] for r in rows],np.int64)
        self.target=F(seed['target'])
        self.goal=seed['target_integer']
        assert self.goal==int(self.target*AREA)+1
        assert int(self.radial_floors.min())>=self.goal

    def read(self,name):
        b=(self.tree/name).read_bytes()
        assert sha(b)==self.manifest[name]
        return b

    def seed_positions(self,part,path):
        with np.load(path,allow_pickle=False) as z:
            pos=z['pos'];score=z['score']
        row=self.seed_rows[part]
        assert pos.dtype==np.int32 and pos.ndim==1 and len(pos)==row[3]
        assert sha(pos.tobytes())==row[-2] and sha(score.tobytes())==row[-1]
        return pos


def run(args):
    work=args.work.resolve()
    work.mkdir(parents=True,exist_ok=True)
    sys.path.insert(0,str(args.computation))
    data=Data(args.computation,args.tree,args.table,args.seed_receipt)
    compare_next=int(F('.83732')*AREA)+1
    C,kernel_sha,binary_sha=kernel(args.computation,work,args.cache)
    parts=sorted(data.seed_rows)
    if args.parts:
        parts=sorted(set(map(int,args.parts.split(','))))
        assert all(p in data.seed_rows for p in parts)
    jobs=sorted(parts,key=lambda p:-data.seed_rows[p][3])
    audits=[]
    for layer,(div,budget) in enumerate(STAGES,1):
        out=work/f'stage{layer}'
        out.mkdir(exist_ok=True)
        def task(part):
            stem=f'part{part:04}'
            path=out/(stem+'.npz')
            receipt=path.with_suffix('.json')
            source=(args.seeds if layer==1 else work/f'stage{layer-1}')/path.name
            if layer==1:
                pos=data.seed_positions(part,source)
            else:
                with np.load(source,allow_pickle=False) as z:
                    pos=z['remaining']
                prior=json.loads(source.with_suffix('.json').read_bytes())
                assert sha(source.read_bytes())==prior['data_sha256']
            assert np.array_equal(pos,np.unique(pos))
            link=sha(source.read_bytes())
            config=dict(part=part,layer=layer,divisor=div,budget=budget,source_sha256=link,
                        kernel_sha256=kernel_sha,goal=data.goal,radial_goal=data.seed['radial_target'])
            if receipt.exists():
                record=json.loads(receipt.read_bytes())
                assert record['config']==config and sha(path.read_bytes())==record['data_sha256']
                return record
            else:
                assert not path.exists()
            with __import__('io').BytesIO(data.read(stem+'.npy')) as b:
                codes=np.load(b,allow_pickle=False)
            assert codes.dtype==np.uint8 and np.all(codes[pos]>=253)
            active=np.zeros(len(codes),np.uint8)
            active[pos]=1
            retained=np.zeros_like(active)
            raw=np.zeros(len(codes),np.int64);clipped=np.zeros_like(raw)
            st=np.zeros(6,np.int64)
            bits=np.array(np.unravel_index(part,(2,)*11),bool)
            lo=np.where(bits,data.mid,data.lo);hi=np.where(bits,data.hi,data.mid)
            tab,counts=data.tables[div]
            assert not C(codes,len(codes),lo,hi,data.hi,active,retained,st,tab,counts,data.goal,48,
                         budget,data.radii,data.radial_floors,raw,clipped)
            assert np.all(retained<=active) and np.all(raw[~active.astype(bool)]==0)
            accepted=pos[raw[pos]>=data.goal]
            pending=np.flatnonzero(retained).astype(np.int32)
            assert np.array_equal(accepted,pos[~retained[pos].astype(bool)])
            assert len(accepted)==st[2] and len(pending)==st[3]
            np.savez_compressed(path,accepted=accepted,remaining=pending,
                                full_floor=raw[accepted],clipped_floor=clipped[accepted])
            record=dict(config=config,data_sha256=sha(path.read_bytes()),counts=st.tolist(),
                        tested=len(pos),accepted=len(accepted),remaining=len(pending),
                        at_next_target=int((raw[accepted]>=compare_next).sum()),
                        capped_by_radial=int(((clipped[accepted]!=0)&(raw[accepted]<clipped[accepted])).sum()))
            receipt.write_text(json.dumps(record,separators=(',',':'))+'\n')
            return record

        st=monotonic()
        values=[]
        with ThreadPoolExecutor(max_workers=args.workers) as pool:
            for i,rr in enumerate(pool.map(task,jobs),1):
                values.append(rr)
                if i%128==0:
                    print(layer,i,len(jobs),'remaining',sum(r['remaining'] for r in values),
                          'seconds',round(monotonic()-st,1),flush=True)
        audit={k:sum(r[k] for r in values) for k in ('tested','accepted','remaining','at_next_target','capped_by_radial')}
        audit['calls']=sum(r['counts'][0] for r in values)
        audit['work_files']=sorted([r['config']['part'],r['data_sha256']] for r in values)
        if layer>1:assert audits[-1]['remaining']==audit['tested']
        audits.append(audit)
        print('stage',layer,{k:v for k,v in audit.items() if k!='work_files'},'seconds',round(monotonic()-st,3),flush=True)
        (work/'stage-ledger.json').write_text(json.dumps(dict(target=str(data.target),integer_goal=data.goal,
                            radial_target=data.seed['radial_target'],parts=len(parts),complete=len(parts)==1024,
                            audits=audits,source_sha=kernel_sha,binary_sha=binary_sha,
                            seed_receipt_sha=sha(args.seed_receipt.read_bytes())),separators=(',',':'))+'\n')
        if layer==args.until:break
    if args.until==4:
        survivors={}
        for part in parts:
            with np.load(work/'stage4'/f'part{part:04}.npz',allow_pickle=False) as z:
                p=z['remaining']
            if len(p):survivors[str(part)]=p.tolist()
        assert sum(map(len,survivors.values()))==audits[-1]['remaining']
        (work/'retained.json').write_text(json.dumps(survivors,separators=(',',':'))+'\n')


if __name__=='__main__':
    parser=ArgumentParser(description=__doc__)
    for n in ('computation','tree','table','seeds','seed_receipt','work','cache'):
        parser.add_argument(n,type=Path)
    parser.add_argument('--workers',type=int,default=4)
    parser.add_argument('--parts')
    parser.add_argument('--until',type=int,choices=(1,2,3,4),default=4)
    args=parser.parse_args()
    assert 1<=args.workers<=12
    run(args)
