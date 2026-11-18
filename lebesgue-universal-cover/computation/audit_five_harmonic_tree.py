"""Intersect one certified 11D exclusion with the original Tomar preorder trees."""
import ast
from fractions import Fraction as F
from hashlib import sha256
from itertools import product
import json
from math import prod


def read_codes(raw):
    assert raw[:6] == b"\x93NUMPY"
    version = tuple(raw[6:8])
    assert version in ((1, 0), (2, 0), (3, 0))
    size = 2 if version == (1, 0) else 4
    start = 8+size
    end = start+int.from_bytes(raw[8:start],"little")
    meta = ast.literal_eval(raw[start:end].decode("utf-8" if version == (3, 0) else "latin1"))
    assert meta["descr"] == "|u1" and meta["fortran_order"] is False
    codes = raw[end:]
    assert meta["shape"] == (len(codes),)
    prefix = [1]
    for code in codes:
        assert prefix[-1]>0 and (code<11 or code>=253)
        prefix.append(prefix[-1]+(1 if code<11 else -1))
    assert prefix[-1]==0
    return codes,prefix


def split(box,axis):
    a,b=box[axis]
    mid=(a+b)/2
    left,right=box.copy(),box.copy()
    left[axis],right[axis]=(a,mid),(mid,b)
    return left,right


def intersection(a,b):
    c=[(max(x,u),min(y,v)) for (x,y),(u,v) in zip(a,b,strict=True)]
    return c if all(x<y for x,y in c) else None


def volume(box):
    return prod(y-x for x,y in box)


def route_child(box,exclusion,point):
    path=[]
    box=box.copy()
    for j in range(11):
        u,v=exclusion[j]
        while box[j][0]<u or box[j][1]>v:
            left,right=split(box,j)
            child=int(point[j]>=right[j][0])
            box=(left,right)[child]
            path.append([j,child])
    return box,path


def remainder_slabs(parent,exclude):
    assert all(p<=a<b<=q for (p,q),(a,b) in zip(parent,exclude,strict=True))
    keep=[]
    for j,(a,b) in enumerate(exclude):
        u,v=parent[j]
        for l,r in ((u,a),(b,v)):
            if l<r:
                keep.append(exclude[:j]+[(l,r)]+parent[j+1:])
    assert volume(parent)==volume(exclude)+sum(map(volume,keep))
    return len(keep)


def audit(tree,center,radius):
    assert len(center)==11 and radius>0
    exclusion=[(x-radius,x+radius) for x in center]
    base=json.loads((tree/"part0000.json").read_bytes())
    assert base["ns"]==[3,5,7,9] and base["split0"]==[2]*11 and base["nparts"]==2048
    root=list(zip(map(F,base["lo"]),map(F,base["hi"]),strict=True))
    assert intersection(root,exclusion)==exclusion
    allowed=[]
    for j,(u,v) in enumerate(root):
        m=(u+v)/2
        allowed.append([(side,a,b) for side,(a,b) in enumerate(((u,m),(m,v)))
                        if max(a,exclusion[j][0])<min(b,exclusion[j][1])])

    cuts=[]
    cut_volume=F(0)
    sources=[]
    candidate_child=None
    for start in product(*allowed):
        part=sum(bit<<(10-j) for j,(bit,_,_) in enumerate(start))
        box=[(a,b) for _,a,b in start]
        hd=tree/f"part{part:04}.json"
        raw_head=hd.read_bytes()
        meta=json.loads(raw_head)
        assert meta["part"]==part
        assert all(meta[k]==base[k] for k in ("ns","split0","nparts","lo","hi"))
        raw_codes=hd.with_suffix(".npy").read_bytes()
        codes,prefix=read_codes(raw_codes)
        assert len(codes)==meta["boxes"]
        sources.append([part,sha256(raw_head).hexdigest(),sha256(raw_codes).hexdigest()])

        queue=[(0,box)]
        while queue:
            index,region=queue.pop()
            cut=intersection(region,exclusion)
            if cut is None:
                continue
            code=codes[index]
            if code<11:
                left,right=split(region,code)
                next_right=prefix.index(prefix[index],index+2)
                queue.extend(((next_right,right),(index+1,left)))
                continue
            removed=volume(cut)
            cut_volume+=removed
            cuts.append([part,index,code,remainder_slabs(region,cut),str(removed/volume(region))])
            if all(a<c<b for (a,b),c in zip(region,center,strict=True)):
                child,path=route_child(region,exclusion,center)
                assert intersection(exclusion,child)==child and candidate_child is None
                candidate_child={
                    "part":part,"pos":index,"extra_splits":path,
                    "ratio_of_original_leaf":str(volume(child)/volume(region)),
                }
    assert cut_volume==volume(exclusion) and candidate_child is not None
    return {
        "source_parts_header_npy_sha256":sources,
        "intersected_original_leaves":len(cuts),
        "certified_box_volume_relative_to_root":str(cut_volume/volume(root)),
        "dyadic_child":candidate_child,
        "cuts_fields":["part","pos","code","retained_slabs","ratio_removed"],
        "cuts":cuts,
    }
