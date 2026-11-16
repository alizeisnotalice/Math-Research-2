#!/usr/bin/env python3
"""Complete true halo/deletion enclosures, exact closures and x-event trees.

Standard library only. Coordinates and masses are exact rational JSON strings.
Every maximum eligible block is evaluated; no sampled radii or receiver window.
"""
from collections import defaultdict
from fractions import Fraction as F
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import sys
import time

PROVENANCE = {
    "closure_method_source":"computation/round064/mass_polytope_variation.py:25-43",
    "source_sha256":"4ef61b815b9c4a9a8972cdaf31716ae7917b59243a29fb8697c8042137788246",
    "use":"coordinate interval mask enumeration and directed integer square-root principle; new self-contained implementation; old code not imported",
    "area_method":"new vertical x-event segment tree, three cover counts; distinct from old horizontal slab merging"
}


def exact(x):
    if isinstance(x,bool)or not isinstance(x,(str,int)):
        raise ValueError("integer or rational string required")
    return F(x)


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def closures(points,masses):
    """All nonempty tight bbox source closures, including ineligible ones."""
    axes=[]
    for axis in (0,1):
        groups=defaultdict(int)
        for i,p in enumerate(points):groups[p[axis]]|=1<<i
        coords=sorted(groups);intervals=[]
        for j in range(len(coords)):
            mask=0
            for k in range(j,len(coords)):
                mask|=groups[coords[k]];intervals.append(mask)
        axes.append(intervals)
    rows=[]
    for mask in sorted({a&b for a in axes[0] for b in axes[1]}-{0}):
        ids=[i for i in range(len(points))if mask>>i&1]
        a=min(points[i][0]for i in ids);b=max(points[i][0]for i in ids)
        c=min(points[i][1]for i in ids);d=max(points[i][1]for i in ids)
        mass=sum(masses[i]for i in ids);span=max(b-a,d-c)
        rows.append({"mask":mask,"ids":ids,"bbox":(a,b,c,d),"mass":mass,
                     "span_squared":span*span,"eligible":mass>span*span})
    return rows


def rectangles(rows,bits):
    """Integer units base*2^bits; exact integer mass units base^2."""
    D=1<<bits;inner=[];outer=[];catalog=[]
    for row in rows:
        if not row["eligible"]:continue
        v=row["mass"]*D*D
        k=math.isqrt(v//4)
        assert 4*k*k<=v<4*(k+1)*(k+1)
        upper=k if 4*k*k==v else k+1
        a,b,c,d=(x*D for x in row["bbox"])
        ri=(b-k,a+k,d-k,c+k)
        ro=(b-upper,a+upper,d-upper,c+upper)
        if ri[0]<ri[1]and ri[2]<ri[3]:inner.append(ri)
        assert ro[0]<ro[1]and ro[2]<ro[3]
        outer.append(ro)
        catalog.append({"mask":row["mask"],"radius_lower_units":k,
                        "radius_upper_units":upper,"inner":ri,"outer":ro})
    return inner,outer,catalog


class TripleCoverTree:
    """Three independent interval-union trees, sharing one y partition."""
    def __init__(self,ys):
        self.ys=ys;self.n=len(ys)-1
        self.tags=[[0]*(4*self.n+8)for _ in range(3)]
        self.lengths=[[0]*(4*self.n+8)for _ in range(3)]

    def change(self,ql,qr,delta,kind,node=1,lo=0,hi=None):
        if hi is None:hi=self.n
        if ql>=hi or qr<=lo:return
        if ql<=lo and hi<=qr:
            for channel in (kind,2):
                self.tags[channel][node]+=delta
                assert self.tags[channel][node]>=0
        else:
            mid=(lo+hi)//2
            self.change(ql,qr,delta,kind,node*2,lo,mid)
            self.change(ql,qr,delta,kind,node*2+1,mid,hi)
        for channel in range(3):
            if self.tags[channel][node]:
                self.lengths[channel][node]=self.ys[hi]-self.ys[lo]
            elif hi-lo==1:self.lengths[channel][node]=0
            else:self.lengths[channel][node]=self.lengths[channel][node*2]+self.lengths[channel][node*2+1]


def sweep_areas(left,right):
    """Return area(left union), area(right union), area(left minus right)."""
    if not left and not right:return (0,0,0)
    ys=sorted({y for rect in left+right for y in rect[2:]})
    yi={y:i for i,y in enumerate(ys)};events=defaultdict(list)
    for kind,rects in enumerate((left,right)):
        for a,b,c,d in rects:
            events[a].append((yi[c],yi[d],1,kind))
            events[b].append((yi[c],yi[d],-1,kind))
    tree=TripleCoverTree(ys);areas=[0,0,0];last=min(events)
    for x,es in sorted(events.items()):
        width=x-last
        for channel in range(3):areas[channel]+=width*tree.lengths[channel][1]
        # At an equal x all events apply before the following open slab.
        for c,d,delta,kind in es:tree.change(c,d,delta,kind)
        last=x
    assert all(tree.lengths[k][1]==0 for k in range(3))
    return areas[0],areas[1],areas[2]-areas[1]


def enc(lo,hi,den):
    assert 0<=lo<=hi
    return {"lower":str(F(lo,den)),"upper":str(F(hi,den)),
            "lower_decimal":float(F(lo,den)),"upper_decimal":float(F(hi,den)),
            "width":str(F(hi-lo,den))}


def evaluate(data,*,bits=72):
    start=time.perf_counter()
    assert bits>=16 and exact(data["threshold"])==1
    if len(data["positions"])!=len(data["masses"]):raise ValueError("array lengths")
    chosen=set(data.get("S_indices",[]))
    pp=[];mm=[];original=[]
    for i,(point,mass)in enumerate(zip(data["positions"],data["masses"])):
        p=tuple(map(exact,point));m=exact(mass)
        if len(p)!=2 or m<0:raise ValueError("invalid atom")
        if m>0:pp.append(p);mm.append(m);original.append(i)
    if not pp or len(set(pp))!=len(pp):raise ValueError("need distinct positive atoms")
    if not chosen.issubset(set(range(len(data["positions"])))):raise ValueError("invalid chosen S")
    S=sum(1<<k for k,i in enumerate(original)if i in chosen)
    base=math.lcm(*(x.denominator for p in pp for x in p),*(m.denominator for m in mm))
    pts=[tuple(int(x*base)for x in p)for p in pp]
    mass=[int(m*base*base)for m in mm]
    full=closures(pts,mass);eligible=[row for row in full if row["eligible"]]
    L=max(row["mass"]for row in eligible)
    maxima=[row for row in eligible if row["mass"]==L]
    # Positivity guarantees a maximum eligible subset is bbox closed: otherwise
    # its closure would remain eligible with strictly larger mass.
    parent=list(range(len(pts)))
    def find(i):
        while parent[i]!=i:parent[i]=parent[parent[i]];i=parent[i]
        return i
    for row in eligible:
        for j in row["ids"][1:]:parent[find(j)]=find(row["ids"][0])
    D=base*(1<<bits);den=D*D
    fi,fo,fc=rectangles(full,bits)
    full_lo=sweep_areas(fi,[])[0];full_hi=sweep_areas(fo,[])[0]
    jobs=[]
    for winner in maxima:
        keep=[i for i in range(len(pts))if not winner["mask"]>>i&1]
        res=closures([pts[i]for i in keep],[mass[i]for i in keep])
        ri,ro,rc=rectangles(res,bits)
        _,rhi,dlo=sweep_areas(fi,ro)
        _,rlo,dhi=sweep_areas(fo,ri)
        own=[i for i in range(len(pts))if winner["mask"]>>i&1]
        ins=closures([pts[i]for i in own],[mass[i]for i in own])
        ii,io,ic=rectangles(ins,bits)
        ilo=sweep_areas(ii,[])[0];ihi=sweep_areas(io,[])[0]
        jobs.append({"maximum_mask":winner["mask"],"maximum_original_indices":[original[i]for i in own],
                     "chosen_S":winner["mask"]==S,"residual_local_to_global":keep,
                     "residual_closure_count":len(res),"residual_eligible_closure_count":sum(r["eligible"]for r in res),
                     "residual_halo":enc(rlo,rhi,den),"internal_halo":enc(ilo,ihi,den),
                     "actual_Delta":enc(dlo,dhi,den),
                     "Delta_over_L":enc(dlo*base*base,dhi*base*base,den*L),
                     "residual_directed_rectangles":rc,"internal_directed_rectangles":ic})
    total=sum(mass);bbox=(min(p[0]for p in pts),max(p[0]for p in pts),min(p[1]for p in pts),max(p[1]for p in pts))
    span=max(bbox[1]-bbox[0],bbox[3]-bbox[2]);gap=F(span*span,base*base)-F(35,33)*F(total,base*base)
    return {"schema":"r072-complete-actual-halo-v1","bits":bits,"coordinate_scale":D,"base_scale":base,
            "positive_source_count":len(pts),"positive_original_indices":original,
            "maximum_mass":str(F(L,base*base)),"maximum_count":len(maxima),"chosen_S_is_maximum":any(row["mask"]==S for row in maxima),
            "full_closure_count":len(full),"full_eligible_closure_count":len(eligible),
            "eligible_component_count":len({find(i)for i in range(len(pts))}),
            "actual_total_mass":str(F(total,base*base)),"actual_span":str(F(span,base)),
            "nonterminal_gap_span_squared_minus_35total_over33":str(gap),"nonterminal":gap>0,
            "full_halo":enc(full_lo,full_hi,den),"all_maximum_block_evaluations":jobs,
            "full_closure_catalog":[{"mask":r["mask"],"bbox_units_base":r["bbox"],"mass_units_base_squared":r["mass"],"span_squared_units_base_squared":r["span_squared"],"eligible":r["eligible"]}for r in full],
            "full_directed_rectangles":fc,"method_provenance":PROVENANCE,
            "wall_seconds":time.perf_counter()-start,
            "scope":"complete finite positive-atom full/residual/internal halos and actual deletion for every maximum eligible block; exact directed enclosures, no window/radius truncation; no general bound"}


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input",type=Path,required=True);ap.add_argument("--output",type=Path,required=True)
    ap.add_argument("--bits",type=int,default=72);args=ap.parse_args()
    assert os.path.ismount('/Volumes/data')
    raw=args.input.read_bytes();out=evaluate(json.loads(raw),bits=args.bits)
    out["receipt"]={"python":sys.version,"platform":platform.platform(),"pid":os.getpid(),
                    "input_sha256":hashlib.sha256(raw).hexdigest(),"script_sha256":sha(Path(__file__))}
    args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:out[k]for k in ("positive_source_count","maximum_count","full_halo","nonterminal","wall_seconds")}))
    print(json.dumps([{k:r[k]for k in ("maximum_mask","chosen_S","actual_Delta")}for r in out["all_maximum_block_evaluations"]]))


if __name__=='__main__':main()
