"""Verify frozen native holdout hosts against the eleven-point test set."""
from argparse import ArgumentParser
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import json

import numpy as np

from certify_multi_host_points import SupportHull, points, q, verify_rows
from tomar_leaf_box import leaf_box

CUT=3684498829194
AREA_BASE=F('.8371')
AREA_TEST=F('.83730')


def native_geometry(root, part, pos, manifest):
    stem=f'part{part:04}'
    js=root/(stem+'.json')
    npy=root/(stem+'.npy')
    assert all(sha256(p.read_bytes()).hexdigest()==manifest[p.name] for p in (js,npy))
    meta=json.loads(js.read_bytes())
    codes=np.load(npy)
    assert meta['part']==part and codes[pos]>=253
    prefix=np.r_[1,1+np.cumsum(np.where(codes<11,1,-1))]
    mid,half=leaf_box(prefix,codes,part,pos,meta)
    f=2**44
    lo=[int((F(a)-F(b))*f) for a,b in zip(mid,half)]
    hi=[int((F(a)+F(b))*f) for a,b in zip(mid,half)]
    assert all(F(i,f)==F(a)-F(b) and F(j,f)==F(a)+F(b) for a,b,i,j in zip(mid,half,lo,hi))
    assert lo[2]<=CUT
    hi[2]=min(hi[2],CUT)
    return lo,hi


def report(tree, source, holdout, point_input):
    manifest_source=(tree/'SHA256SUMS').read_bytes()
    assert sha256(manifest_source).hexdigest()==holdout['manifest_sha256']
    assert sha256(source).hexdigest()==holdout['source_sha256']
    assert sha256(point_input).hexdigest()==holdout['point_certificate_sha256']
    manifest={f.removeprefix('./'):h for h,f in (l.split() for l in manifest_source.decode().splitlines())}
    leaf_set=json.loads(source)
    cert=json.loads(point_input)
    P,_=points(cert)
    d=holdout['denominator']
    assert d==2**44 and holdout['seed']==20261118
    sample=np.random.default_rng(holdout['seed'])
    selected=[]
    for key,ids in sorted(leaf_set.items(),key=lambda d:int(d[0])):
        part=int(key)
        vals=list(map(int,sample.choice(ids,min(len(ids),max(1,32//len(leaf_set))),replace=False)))
        if part==180 and 23955 not in vals:
            vals.append(23955)
        selected.extend((part,i) for i in sorted(vals))
    assert selected==[(x['part'],x['leaf']) for x in holdout['rows']]
    full={}
    centre={}
    for row in holdout['rows']:
        part,leaf=row['part'],row['leaf']
        lo,hi=native_geometry(tree,part,leaf,manifest)
        assert lo==row['lo'] and hi==row['hi']
        a,b=([F(x,d) for x in z] for z in (lo,hi))
        nom=[(p+t)/2 for p,t in zip(a,b)]
        hw=[(t-p)/2 for p,t in zip(a,b)]
        for tag,spread,out in [('centre',None,centre),('full',hw,full)]:
            if row[tag] is None:
                continue
            target=F(row[tag]['target'])
            assert target>0
            host=SupportHull([q(x) for x in nom],spread)
            views={}
            for sign in (1,-1):
                views[str(sign)]=verify_rows(host,P,sign,row[tag]['hands'][str(sign)],target)
            result={'padding':str(target),'hands':views}
            if tag=='full':
                delta=target**2/(4*AREA_TEST)
                assert AREA_BASE+delta<AREA_TEST
                result['local_six_input_area_lower']=str(AREA_BASE+delta)
            out[f'{part}/{leaf}']=result
    return {'sample_size':len(selected),'strict_midpoint_obstructions':len(centre),
            'strict_whole_native_leaf_obstructions':len(full),'centre':centre,'full':full}


def main():
    ap=ArgumentParser(description=__doc__)
    ap.add_argument('tree',type=Path)
    ap.add_argument('remaining_leaf_snapshot',type=Path)
    ap.add_argument('holdout_data',type=Path)
    ap.add_argument('point_certificate',type=Path)
    args=ap.parse_args()
    data=args.holdout_data.read_bytes()
    out=report(args.tree,args.remaining_leaf_snapshot.read_bytes(),json.loads(data),args.point_certificate.read_bytes())
    out['holdout_data_sha256']=sha256(data).hexdigest()
    out['script_sha256']=sha256(Path(__file__).read_bytes()).hexdigest()
    print(json.dumps(out,indent=2))


if __name__=='__main__':
    main()
