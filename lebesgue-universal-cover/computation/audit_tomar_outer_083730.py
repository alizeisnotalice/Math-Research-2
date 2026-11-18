"""Replay score reuse, radial cuts and rigorous deficits for the 0.83730 outer domain."""
from argparse import ArgumentParser
from fractions import Fraction as F
from hashlib import sha256
from io import BytesIO
from math import isqrt
from pathlib import Path
from zipfile import ZipFile
import ctypes
import json
import os
import subprocess
import sys

import numpy as np

FINE=2**44
AREA=2**57
TARGET=F('.83730')
OLD=F('.837295006222')
RADII=tuple(map(F,('.06835','.06650','.06004','.05699')))
SMALL_DEFICITS=tuple(map(F,('.0000001','.0000005','.000001','.000002')))
FILES={'score_index':'tomar-grouped-08371-index.json','score_ledger':'tomar-grouped-08371-ledger.json',
       'slice':'five-orbit-slice.json','localization':'five-global-localization.json',
       'radial':'five-body-radial-affine.json','kernel':'tomar_grouped_scored_fan.cpp'}

CPP=r"""
#include <algorithm>
#include <array>
#include <cstdint>
#include <vector>
using L=std::int64_t;
using W=__int128_t;
using Box=std::array<L,11>;
struct Filter {
    const std::uint8_t* codes;
    const std::vector<int>& stop;
    const std::vector<int>& prefix;
    const std::uint8_t* selected;
    const L* r;
    std::uint8_t* output;
    L excluded=0, remaining=0;
    bool scan(int at,Box& lo,Box& hi,int depth) {
        if(depth>256)return false;
        int n=prefix[stop[at]]-prefix[at];
        if(!n)return true;
        const auto distance=[&](int j){return lo[j]>0?lo[j]:hi[j]<0?-hi[j]:0;};
        for(int b=0;b<4;++b){
            const L x=distance(3*b), y=distance(3*b+1);
            if(W(x)*x+W(y)*y>W(r[b])*r[b]){excluded+=n;return true;}
        }
        const int k=codes[at];
        if(k>=253){
            if(!selected[at])return false;
            output[at]=1;++remaining;return true;
        }
        const L a=lo[k],b=hi[k];
        if((a+b)%2)return false;
        const L m=(a+b)/2;
        hi[k]=m;
        if(!scan(at+1,lo,hi,depth+1))return false;
        hi[k]=b;lo[k]=m;
        if(!scan(stop[at+1],lo,hi,depth+1))return false;
        lo[k]=a;
        return true;
    }
};
extern "C" int audit_outer_radial(const std::uint8_t* codes,int n,const L* low,const L* high,
                                  const L* r,const std::uint8_t* subset,std::uint8_t* keep,L* counts) {
    if(n<1 || n>0x10000000)return 1;
    std::vector<int> stop(n),prefix(n+1);
    for(int i=n-1;i>=0;--i) {
        if(codes[i]>=253)stop[i]=i+1;
        else {
            if(codes[i]>=11 || i+1>=n || stop[i+1]>=n)return 2;
            stop[i]=stop[stop[i+1]];
        }
    }
    if(stop[0]!=n)return 3;
    for(int i=0;i<n;++i){
        if(subset[i]>1 || (subset[i] && codes[i]<253))return 4;
        prefix[i+1]=prefix[i]+subset[i];
    }
    Box a,b; std::copy(low,low+11,a.begin()); std::copy(high,high+11,b.begin());
    Filter f{codes,stop,prefix,subset,r,keep};
    if(!f.scan(0,a,b,0) || f.remaining+f.excluded!=prefix[n])return 5;
    counts[0]=f.excluded;counts[1]=f.remaining;
    return 0;
}
"""


def sha(data):
    return sha256(data).hexdigest()


def ceil(x):
    return -(-x.numerator//x.denominator)


def score_floor(b):
    assert 0<=b<=65535
    quantized=F(8371000+b-1,10**7) if b else F('.8371')
    return ceil(AREA*quantized)


def support_surplus(amount):
    """One safe outward support surplus on the 10^-8 grid for a deficit amount."""
    need=4*TARGET*amount*10**16
    candidate=isqrt(ceil(need))
    if candidate*candidate<need:
        candidate+=1
    assert F(candidate,10**8)**2>=4*TARGET*amount
    return F(candidate,10**8)


def compile_filter(work,cache):
    source=work/'outer-radial.cpp'
    lib=work/'outer-radial.dylib'
    source.write_text(CPP)
    compiler=os.environ.get('CXX','clang++')
    subprocess.run([compiler,'-O3','-std=c++17','-shared','-fPIC','-ftrapv',
                    f'-fmodules-cache-path={cache}',str(source),'-o',str(lib)],check=True)
    C=ctypes.CDLL(str(lib))
    v=lambda t:np.ctypeslib.ndpointer(dtype=t,ndim=1,flags='C_CONTIGUOUS')
    C.audit_outer_radial.argtypes=[v(np.uint8),ctypes.c_int,*[v(np.int64)]*3,
                                  v(np.uint8),v(np.uint8),v(np.int64)]
    return C.audit_outer_radial


def run(args):
    sys.path.insert(0,str(args.computation))
    from tomar_d3_mask import leaf_mask, CUT
    from flint import arb,ctx

    ctx.prec=160
    raw={k:(args.computation/v).read_bytes() for k,v in FILES.items()}
    docs={k:json.loads(v) for k,v in raw.items() if k!='kernel'}
    idx,ledger,S,loc,radial=(docs[k] for k in ('score_index','score_ledger','slice','localization','radial'))
    assert idx['audit_sha256']==sha(raw['score_ledger'])
    assert ledger['slice_index_sha256']==sha(raw['slice'])
    assert idx['complete'] and ledger['complete'] and idx['target']==ledger['target']=='0.8371'
    assert len(idx['archives'])==4 and not ledger['historical_area_inputs']
    assert loc['denominator']==FINE and F(loc['target'])==OLD
    assert len(loc['retained'])==19 and sum(map(len,loc['retained'].values()))==501

    manifest=(args.tree/'SHA256SUMS').read_bytes()
    assert sha(manifest)==S['manifest']==ledger['input_tree_sha256']
    files={n.removeprefix('./'):v for v,n in (line.split() for line in manifest.decode().splitlines())}
    rootname='part0000.json'
    root_raw=(args.tree/rootname).read_bytes()
    assert sha(root_raw)==files[rootname]
    root=json.loads(root_raw)
    assert root['ns']==[3,5,7,9] and root['split0']==[2]*11
    minbox=[F(x)*FINE for x in root['lo']]
    maxbox=[F(x)*FINE for x in root['hi']]
    assert all(x.denominator==1 for x in minbox+maxbox)
    LL=np.array(list(map(int,minbox)),np.int64);HH=np.array(list(map(int,maxbox)),np.int64)
    assert np.all((LL+HH)%2==0)
    mitt=(LL+HH)//2
    assert arb.pi()/15 < arb(CUT.numerator)/CUT.denominator and F(root['hi'][2])>=CUT
    for n,k in ((7,5),(9,8)):
        z=F(root['hi'][k])
        assert arb(z.numerator)/z.denominator > 2*arb.pi()/n
    assert all(F(root['lo'][j])<0<F(root['hi'][j]) for j in (0,1,3,4,6,7,9,10))
    assert all(arb(x)>arb('.83730') for x in ledger['root_exterior'])

    assert (radial['source_sha256']==sha((args.computation/'certify_radial_affine_fans.py').read_bytes())
            and radial['table_sha']==sha((args.computation/'five-body-radial-material.npz').read_bytes()))
    radii=np.array([ceil(r*FINE) for r in RADII],np.int64)
    floors=[]
    for j,r in enumerate(radial['rows']):
        assert F(r['radius'])==RADII[j]
        bound=F(r['uniform_twice_area_numerator'],AREA)
        assert bound>TARGET
        floors.append(bound)
        for coord in (3*j,3*j+1):
            assert minbox[coord]<-radii[j]<radii[j]<maxbox[coord]
    if args.replay_radial:
        result=subprocess.check_output([sys.executable,'-I','-B',args.computation/'certify_radial_affine_fans.py'])
        assert json.loads(result)==radial

    row_by_part={p:row for p,*row in ledger['receipt_hashes']}
    assert len(row_by_part)==1024
    args.work.mkdir(parents=True,exist_ok=True)
    cast=compile_filter(args.work,args.cache)
    if args.emit is not None:
        args.emit.mkdir(parents=True,exist_ok=True)
    G,G_0=(ceil(v*AREA) for v in (TARGET,OLD))
    assert F(G,AREA)>TARGET
    low=np.array([score_floor(b) for b in range(65536)],dtype=np.int64)
    wanted=int(np.searchsorted(low,G))
    assert wanted==2001
    histogram=np.zeros(65536,dtype=np.int64)
    limits=(*SMALL_DEFICITS,F(G-G_0,AREA))
    ranges=np.zeros(len(limits),dtype=np.int64)
    thresholds=np.array([int(v*AREA) for v in limits],dtype=np.int64)
    outer_stronger=0
    D=np.zeros(7,dtype=np.int64)
    done=set()
    results=[]
    for archive in idx['archives']:
        ar=args.computation/archive['name']
        assert len(ar.read_bytes())==archive['bytes'] and sha(ar.read_bytes())==archive['sha256']
        with ZipFile(ar) as z:
            for member in z.namelist():
                if not (member.startswith('native/') and member.endswith('.npz')):
                    continue
                stem=Path(member).stem
                part=int(stem[4:])
                assert part not in done
                done.add(part)
                packed=z.read(member); receipt=z.read(member.replace('.npz','.json'))
                assert [sha(receipt),sha(packed)]==row_by_part[part][:2]
                native=json.loads(receipt)
                source=(args.tree/(stem+'.npy')).read_bytes()
                head=(args.tree/(stem+'.json')).read_bytes()
                assert sha(source)==files[stem+'.npy']==native['codes_sha256']
                assert sha(head)==files[stem+'.json']==native['header_sha256']
                metadata=json.loads(head)
                assert all(metadata[k]==root[k] for k in ('lo','hi','split0','ns'))
                codes=np.load(BytesIO(source),allow_pickle=False)
                cls=(leaf_mask(codes,part,F(root['hi'][2])) if args.masks is None
                     else np.load(args.masks/(stem+'.npy'),allow_pickle=False))
                assert cls.dtype==np.uint8 and len(cls)==len(codes)
                assert sha(cls.tobytes())==S['rows'][part][-1]==native['slice_mask_sha256']
                pos=np.flatnonzero((cls==1)|(cls==2)).astype(np.int32)
                with np.load(BytesIO(packed),allow_pickle=False) as a:
                    score=a['scores']
                assert score.dtype==np.uint16 and score.shape==pos.shape and len(pos)==native['all_original']
                histogram+=np.bincount(score,minlength=65536)
                test=score<wanted
                request=pos[test]
                active=np.zeros(len(codes),dtype=np.uint8)
                active[request]=1
                saved=np.zeros_like(active)
                stats=np.zeros(2,np.int64)
                midbits=np.array(np.unravel_index(part,(2,)*11),bool)
                l=np.where(midbits,mitt,LL)
                h=np.where(midbits,HH,mitt)
                assert not cast(codes,len(codes),l,h,radii,active,saved,stats)
                kept=np.flatnonzero(saved).astype(np.int32)
                assert np.all(saved<=active) and stats.sum()==len(request) and len(kept)==stats[1]
                ss=np.searchsorted(pos,kept)
                assert np.array_equal(pos[ss],kept)
                if args.compare_radial is not None:
                    past=np.load(args.compare_radial/(stem+'.npz'),allow_pickle=False)['pos']
                    assert np.array_equal(kept,request[np.isin(request,past,assume_unique=True)])
                high=np.asarray(loc['retained'].get(str(part),[]),dtype=np.int32)
                is_inner=np.isin(kept,high,assume_unique=True)
                assert np.array_equal(kept[is_inner],high)
                inherited=np.maximum(G_0,low[score[ss[~is_inner]]])
                outer_stronger+=int((inherited>G_0).sum())
                deficit=G-inherited
                assert np.all(deficit>0)
                tier=np.searchsorted(thresholds,deficit,side='left')
                assert np.all(tier<len(limits))
                ranges+=np.bincount(tier,minlength=len(limits))
                group=[len(pos),int((~test).sum()),int(stats[0]),len(kept),len(high),len(deficit),
                       int((score[ss]==0).sum())]
                D+=group
                results.append([part,*group,sha(kept.tobytes()),sha(score[ss].tobytes())])
                if args.emit is not None:
                    np.savez_compressed(args.emit/(stem+'.npz'),pos=kept,
                                        score=score[ss],inside501=is_inner)
    assert done==set(row_by_part) and int(D[0])==6710314 and D[1]+D[2]+D[3]==D[0]
    assert D[4]==501 and D[4]+D[5]==D[3] and int(ranges.sum())==D[5]
    print(json.dumps(dict(target=str(TARGET),old_target=str(OLD),int_target=G,int_old_target=G_0,
        smallest_score_accept_margin=str(F(G,AREA)-TARGET),
        inherited_outer_loss_upper=str(F(G-G_0,AREA)),
        minimum_radial_margin=str(min(floors)-TARGET),
        qualifying_score=wanted,exact_target_bucket=int(histogram[wanted]),
        totals=dict(zip(('canonical','score_accept','radial_accept_after_score','remaining',
                         'within_old_501','remaining_outer','unscored_remaining'),map(int,D))),
        outer_with_stronger_score=outer_stronger,outer_using_old_threshold_only=int(D[5])-outer_stronger,
        outer_deficit_bins=[dict(max_deficit=str(b),leaf_count=int(n),
                     protrusion_sufficient=str(support_surplus(b))) for b,n in zip(limits,ranges)],
        provenance={k:sha(v) for k,v in raw.items()},tree_manifest_sha256=sha(manifest),
        script_sha256=sha(Path(__file__).read_bytes()),rows=sorted(results)),separators=(',',':')))


if __name__=='__main__':
    ap=ArgumentParser(description=__doc__)
    for n in ('computation','tree','work','cache'):
        ap.add_argument(n,type=Path)
    ap.add_argument('--masks',type=Path,help='reuse exact masks with verified per-part digests')
    ap.add_argument('--emit',type=Path,help='write per-part exact survivor arrays')
    ap.add_argument('--compare-radial',type=Path,help='independent radial snapshot')
    ap.add_argument('--replay-radial',action='store_true')
    a=ap.parse_args()
    run(a)
