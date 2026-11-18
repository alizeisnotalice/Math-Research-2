"""Exact original five-body mask for a new target, reusing authenticated old area integers."""
from argparse import ArgumentParser
from fractions import Fraction as F
from io import BytesIO
from pathlib import Path
from zipfile import ZipFile
import json
import sys

import numpy as np
from flint import arb,ctx

from audit_tomar_outer_083730 import AREA, FINE, ceil, compile_filter, score_floor, sha


def q(x):
    f=F(x)
    return arb(f.numerator)/f.denominator


def radius_rows(prior,required):
    output=[]
    for row in prior['rows']:
        n=row['n'];r=F(row['radius']);value=F(row['uniform_twice_area_numerator'],AREA)
        alpha=arb.pi()/(2*n)
        base=arb.pi()/4+n*(alpha.tan()-alpha)/4
        bound=F(int(float(base.mid())*10**10)+1,10**10)
        assert base<q(bound) and bound<value
        if value>required:
            outer=ceil(r*FINE)
        else:
            outer=ceil(r*FINE*(required-bound)/(value-bound))+1
        new=value+(F(outer,FINE)/r-1)*(value-bound)
        lower=int(new*AREA)
        assert F(outer,FINE)>=r and F(lower,AREA)>required
        output.append({'n':n,'previous_radius':str(r),'old_floor':row['uniform_twice_area_numerator'],
                       'center_upper':str(bound),'integer_radius':outer,'integer_area_lower':lower})
    assert [x['n'] for x in output]==[3,5,7,9]
    return output


def select(args):
    ctx.prec=160
    main=args.computation
    sys.path.insert(0,str(main))
    from tomar_d3_mask import leaf_mask, CUT

    paths={k:main/f for k,f in {
        'score_index':'tomar-grouped-08371-index.json','score_ledger':'tomar-grouped-08371-ledger.json',
        'slice_index':'five-orbit-slice.json','radial_certificate':'five-body-radial-affine.json'}.items()}
    raw={k:p.read_bytes() for k,p in paths.items()}
    docs={k:json.loads(v) for k,v in raw.items()}
    idx=docs['score_index'];ledger=docs['score_ledger'];S=docs['slice_index'];rad=docs['radial_certificate']
    assert sha(raw['score_ledger'])==idx['audit_sha256'] and sha(raw['slice_index'])==ledger['slice_index_sha256']
    assert idx['target']==ledger['target']=='0.8371' and idx['complete'] and ledger['complete']
    assert not ledger['historical_area_inputs']
    assert rad['source_sha256']==sha((main/'certify_radial_affine_fans.py').read_bytes())
    assert rad['table_sha']==sha((main/'five-body-radial-material.npz').read_bytes())
    assert args.radius_target>=args.target>F('.8371')
    rr=radius_rows(rad,args.radius_target)
    rr_int=np.array([row['integer_radius'] for row in rr],np.int64)
    goal=int(args.target*AREA)+1
    wanted=int(np.searchsorted(np.array([score_floor(b) for b in range(65536)],np.int64),goal))
    assert 0<wanted<65536 and F(goal,AREA)>args.target
    manifest=(args.tree/'SHA256SUMS').read_bytes()
    assert sha(manifest)==ledger['input_tree_sha256']==S['manifest']
    names={n.removeprefix('./'):h for h,n in (s.split() for s in manifest.decode().splitlines())}
    base=(args.tree/'part0000.json').read_bytes()
    assert sha(base)==names['part0000.json']
    root=json.loads(base)
    lo,hi=(np.array([int(F(v)*FINE) for v in root[k]],np.int64) for k in ('lo','hi'))
    assert all(F(i,FINE)==F(v) for k,a in [('lo',lo),('hi',hi)] for i,v in zip(a,root[k],strict=True))
    assert np.all((lo+hi)%2==0)
    mid=(lo+hi)//2
    for i,R in enumerate(rr_int):
        for j in (3*i,3*i+1):
            assert lo[j]<-R<R<hi[j]
    assert arb.pi()/15<q(CUT)<=q(F(int(hi[2]),FINE))
    assert all(q(F(int(hi[i]),FINE))>2*arb.pi()/n for i,n in ((5,7),(8,9)))
    assert all(lo[i]==0 for i in (2,5,8)) and root['split0']==[2]*11
    assert all(arb(q0)>q(args.target) for q0 in ledger['root_exterior'])
    args.work.mkdir(parents=True,exist_ok=True)
    radial_filter=compile_filter(args.work,args.cache)
    if args.emit is not None:
        args.emit.mkdir(parents=True,exist_ok=True)
    rows={i:x for i,*x in ledger['receipt_hashes']}
    seen=set();records=[];stats=np.zeros(5,np.int64)
    for a in idx['archives']:
        path=main/a['name']
        assert sha(path.read_bytes())==a['sha256']
        with ZipFile(path) as z:
            for n in z.namelist():
                if not (n.startswith('native/') and n.endswith('.npz')):
                    continue
                stem=Path(n).stem;part=int(stem[4:])
                assert part not in seen
                seen.add(part)
                pkt,receipt=z.read(n),z.read(n.replace('.npz','.json'))
                assert [sha(receipt),sha(pkt)]==rows[part][:2]
                meta=json.loads(receipt)
                orig,head=((args.tree/(stem+suffix)).read_bytes() for suffix in ('.npy','.json'))
                assert sha(orig)==names[stem+'.npy']==meta['codes_sha256']
                assert sha(head)==names[stem+'.json']==meta['header_sha256']
                assert all(json.loads(head)[k]==root[k] for k in ('lo','hi','ns','split0'))
                codes=np.load(BytesIO(orig),allow_pickle=False)
                labels=(leaf_mask(codes,part,F(root['hi'][2])) if args.masks is None
                        else np.load(args.masks/(stem+'.npy'),allow_pickle=False))
                assert labels.dtype==np.uint8 and len(labels)==len(codes)
                assert sha(labels.tobytes())==S['rows'][part][-1]==meta['slice_mask_sha256']
                pos=np.flatnonzero((labels==1)|(labels==2)).astype(np.int32)
                with np.load(BytesIO(pkt),allow_pickle=False) as data:
                    score=data['scores']
                assert score.dtype==np.uint16 and len(score)==len(pos)==meta['all_original']
                weak=score<wanted
                rawpos=pos[weak]
                active=np.zeros(len(codes),np.uint8);active[rawpos]=1
                keep=np.zeros_like(active)
                count=np.zeros(2,np.int64)
                bits=np.array(np.unravel_index(part,(2,)*11),bool)
                lo_root=np.where(bits,mid,lo);hi_root=np.where(bits,hi,mid)
                assert not radial_filter(codes,len(codes),lo_root,hi_root,rr_int,active,keep,count)
                assert count.sum()==len(rawpos) and np.all(keep<=active)
                out=np.flatnonzero(keep).astype(np.int32)
                assert len(out)==count[1]
                zscore=score[np.searchsorted(pos,out)]
                counts=[len(pos),int((~weak).sum()),int(count[0]),len(out),int((zscore==0).sum())]
                stats+=counts
                records.append([part,*counts,sha(out.tobytes()),sha(zscore.tobytes())])
                if args.emit is not None:
                    np.savez_compressed(args.emit/(stem+'.npz'),pos=out,score=zscore)
    assert seen==set(rows) and len(seen)==1024 and stats[0]==6710314
    assert stats[1]+stats[2]+stats[3]==stats[0]
    return dict(target=str(args.target),target_integer=goal,radial_target=str(args.radius_target),
                min_bucket=wanted,radii=rr,totals=dict(zip(
                    ('canonical','score_cover','radial_extra','remaining','unscored_remaining'),map(int,stats))),
                provenance={k:sha(v) for k,v in raw.items()},tree_manifest_sha256=sha(manifest),
                rows=sorted(records),script_sha256=sha(Path(__file__).read_bytes()))


if __name__=='__main__':
    p=ArgumentParser(description=__doc__)
    for name in ('computation','tree','work','cache'):
        p.add_argument(name,type=Path)
    p.add_argument('--target',type=F,default=F('.83731'))
    p.add_argument('--radius-target',type=F,default=F('.83732'))
    p.add_argument('--masks',type=Path)
    p.add_argument('--emit',type=Path)
    print(json.dumps(select(p.parse_args()),separators=(',',':')))
