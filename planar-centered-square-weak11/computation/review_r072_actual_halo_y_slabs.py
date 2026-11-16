#!/usr/bin/env python3
"""Independent complete halo/deletion audit from atoms, direct y-slab merging."""
import argparse
from fractions import Fraction as F
from itertools import combinations_with_replacement, product
import hashlib
import json
from math import isqrt, lcm
from pathlib import Path
import sys
import time


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def direct_closures(points,masses):
    """Direct geometric coordinate rectangles and literal membership scans.

    No interval masks are constructed. Canonical tight bboxes deduplicate
    candidate rectangles whose occupied coordinate extremes are smaller.
    """
    if not points:return [],0
    xs=sorted({p[0] for p in points});ys=sorted({p[1] for p in points})
    found={};scans=0
    for (a,b),(c,d) in product(combinations_with_replacement(xs,2),combinations_with_replacement(ys,2)):
        scans+=1
        ids=tuple(i for i,(x,y) in enumerate(points) if a<=x<=b and c<=y<=d)
        if not ids or ids in found:continue
        bb=(min(points[i][0] for i in ids),max(points[i][0] for i in ids),min(points[i][1] for i in ids),max(points[i][1] for i in ids))
        assert ids==tuple(i for i,(x,y) in enumerate(points) if bb[0]<=x<=bb[1] and bb[2]<=y<=bb[3])
        mass=sum((masses[i] for i in ids),F(0));diam=max(bb[1]-bb[0],bb[3]-bb[2])
        found[ids]=dict(mask=sum(1<<i for i in ids),ids=ids,bbox=bb,mass=mass,eligible=diam**2<mass,span_squared=diam**2)
    return sorted(found.values(),key=lambda r:r['mask']),scans


def receivers(rows,scale):
    inner=[];outer=[];catalog=[]
    for row in rows:
        if not row['eligible']:continue
        mass=row['mass'];sq=mass*scale**2/4
        lower=isqrt(sq.numerator//sq.denominator)
        assert F(lower*lower)<=sq<F((lower+1)**2)
        upper=lower if lower*lower==sq else lower+1
        a,b,c,d=[p*scale for p in row['bbox']];assert all(p.denominator==1 for p in (a,b,c,d))
        a,b,c,d=map(int,(a,b,c,d))
        ri=(b-lower,a+lower,d-lower,c+lower);ro=(b-upper,a+upper,d-upper,c+upper)
        if ri[0]<ri[1] and ri[2]<ri[3]:inner.append(ri)
        assert ro[0]<ro[1] and ro[2]<ro[3];outer.append(ro)
        catalog.append(dict(mask=row['mask'],radius_lower_units=lower,radius_upper_units=upper,inner=list(ri),outer=list(ro)))
    return inner,outer,catalog


def merge(intervals):
    result=[]
    for a,b in sorted(intervals):
        if result and a<=result[-1][1]:result[-1]=(result[-1][0],max(result[-1][1],b))
        else:result.append((a,b))
    return result


def difference_length(left,right):
    # Subtract disjoint merged right intervals from disjoint left intervals.
    total=0;j=0
    for a,b in left:
        cursor=a
        while j<len(right) and right[j][1]<=a:j+=1
        k=j
        while k<len(right) and right[k][0]<b:
            lo,hi=right[k]
            if cursor<lo:total+=min(b,lo)-cursor
            cursor=max(cursor,min(b,hi))
            if cursor==b:break
            k+=1
        total+=b-cursor
    return total


def y_slabs(left,right):
    """Literal horizontal slab scans; no sweep tree and no producer code."""
    endpoints=sorted({y for r in left+right for y in r[2:]})
    fa=ra=delta=0;slabs=0;rectangle_membership_tests=0
    for lo,hi in zip(endpoints,endpoints[1:]):
        mid2=lo+hi
        li=merge([(a,b) for a,b,c,d in left if 2*c<mid2<2*d])
        ri=merge([(a,b) for a,b,c,d in right if 2*c<mid2<2*d])
        fa+=(hi-lo)*sum(b-a for a,b in li);ra+=(hi-lo)*sum(b-a for a,b in ri)
        delta+=(hi-lo)*difference_length(li,ri)
        slabs+=1;rectangle_membership_tests+=len(left)+len(right)
    return (fa,ra,delta),dict(y_slabs=slabs,literal_rectangle_membership_tests=rectangle_membership_tests)


def compare_interval(saved,lo,hi,scale):
    assert F(saved['lower'])==F(lo,scale**2) and F(saved['upper'])==F(hi,scale**2)
    return dict(lower=str(F(lo,scale**2)),upper=str(F(hi,scale**2)),display_only=float(F(lo+hi,2*scale**2)))


def main():
    ap=argparse.ArgumentParser()
    for name in ('input','saved','producer','output'):ap.add_argument('--'+name,type=Path,required=True)
    ap.add_argument('--bits',type=int,default=72);a=ap.parse_args();start=time.perf_counter()
    data=json.loads(a.input.read_text());saved=json.loads(a.saved.read_text())
    assert F(data['threshold'])==1 and saved['bits']==a.bits
    assert saved['receipt']['input_sha256']==sha(a.input) and saved['receipt']['script_sha256']==sha(a.producer)
    positions=[];masses=[];original=[]
    for i,(q,m) in enumerate(zip(data['positions'],data['masses'])):
        assert F(m)>=0
        if F(m)>0:positions.append(tuple(map(F,q)));masses.append(F(m));original.append(i)
    assert len(set(positions))==len(positions) and positions
    assert len(positions)==saved['positive_source_count'] and original==saved['positive_original_indices']
    base=lcm(*(c.denominator for q in positions for c in q),*(m.denominator for m in masses));scale=base*(1<<a.bits)
    assert base==saved['base_scale'] and scale==saved['coordinate_scale']
    full,scans=direct_closures(positions,masses);eligible=[r for r in full if r['eligible']];L=max(r['mass'] for r in eligible)
    winners=[r for r in eligible if r['mass']==L];assert L==F(saved['maximum_mass'])
    assert len(full)==saved['full_closure_count'] and len(eligible)==saved['full_eligible_closure_count']
    actual_catalog=[dict(mask=r['mask'],bbox_units_base=[int(x*base) for x in r['bbox']],mass_units_base_squared=int(r['mass']*base**2),span_squared_units_base_squared=int(r['span_squared']*base**2),eligible=r['eligible']) for r in full]
    assert actual_catalog==saved['full_closure_catalog']
    adjacency=[set() for _ in positions]
    for row in eligible:
        for i in row['ids'][1:]:adjacency[row['ids'][0]].add(i);adjacency[i].add(row['ids'][0])
    unseen=set(range(len(positions)));components=0
    while unseen:
        components+=1;stack=[unseen.pop()]
        while stack:
            fresh=adjacency[stack.pop()]&unseen;unseen-=fresh;stack.extend(fresh)
    assert components==saved['eligible_component_count']
    fi,fo,fc=receivers(full,scale);assert fc==saved['full_directed_rectangles']
    (fl,_,_),st1=y_slabs(fi,[]);(fh,_,_),st2=y_slabs(fo,[]);full_area=compare_interval(saved['full_halo'],fl,fh,scale)
    jobmap={r['maximum_mask']:r for r in saved['all_maximum_block_evaluations']};assert set(jobmap)=={r['mask'] for r in winners}
    chosen={original.index(i) for i in data.get('S_indices',[]) if i in original};jobs=[];all_stats=[st1,st2]
    for w in winners:
        job=jobmap[w['mask']];keep=[i for i in range(len(positions)) if i not in w['ids']]
        assert keep==job['residual_local_to_global'] and [original[i] for i in w['ids']]==job['maximum_original_indices']
        assert job['chosen_S']==(set(w['ids'])==chosen)
        residual,rs=direct_closures([positions[i] for i in keep],[masses[i] for i in keep]);internal,ins=direct_closures([positions[i] for i in w['ids']],[masses[i] for i in w['ids']])
        ri,ro,rc=receivers(residual,scale);ii,io,ic=receivers(internal,scale)
        assert len(residual)==job['residual_closure_count'] and sum(r['eligible'] for r in residual)==job['residual_eligible_closure_count']
        assert rc==job['residual_directed_rectangles'] and ic==job['internal_directed_rectangles']
        (_,rhi,dlo),s1=y_slabs(fi,ro);(_,rlo,dhi),s2=y_slabs(fo,ri)
        (ilo,_,_),s3=y_slabs(ii,[]);(ihi,_,_),s4=y_slabs(io,[]);all_stats +=[s1,s2,s3,s4]
        delta=compare_interval(job['actual_Delta'],dlo,dhi,scale);rh=compare_interval(job['residual_halo'],rlo,rhi,scale);ih=compare_interval(job['internal_halo'],ilo,ihi,scale)
        assert F(job['Delta_over_L']['lower'])==F(dlo,scale**2)/L and F(job['Delta_over_L']['upper'])==F(dhi,scale**2)/L
        jobs.append(dict(maximum_mask=w['mask'],chosen_S=job['chosen_S'],mass=str(L),actual_Delta=delta,residual_halo=rh,internal_halo=ih,residual_direct_coordinate_rectangles=rs,internal_direct_coordinate_rectangles=ins))
    span=max(max(q[k] for q in positions)-min(q[k] for q in positions) for k in (0,1));gap=span**2-F(35,33)*sum(masses,F(0))
    assert span==F(saved['actual_span']) and sum(masses,F(0))==F(saved['actual_total_mass'])
    assert gap==F(saved['nonterminal_gap_span_squared_minus_35total_over33']) and saved['nonterminal']==(gap>0)
    result=dict(status='PASS_INDEPENDENT_FULL_ATOMS_HORIZONTAL_SLAB_AREA_AUDIT',python=sys.version,source_sha256=sha(Path(__file__)),input_sha256=sha(a.input),saved_result_sha256=sha(a.saved),producer_sha256=sha(a.producer),bits=a.bits,base_scale=base,positive_sources=len(positions),direct_full_coordinate_rectangles=scans,full_closures=len(full),full_eligible_closures=len(eligible),maximum_mass=str(L),maximum_count=len(winners),full_halo=full_area,all_maximum_blocks=jobs,nonterminal_gap=str(gap),actual_slabs=sum(s['y_slabs'] for s in all_stats),actual_literal_rectangle_membership_tests=sum(s['literal_rectangle_membership_tests'] for s in all_stats),all_exact_catalog_and_area_endpoints_matched=True,elapsed_seconds=time.perf_counter()-start,scope='Actual raw atoms, direct coordinate rectangles and literal membership dedup, all eligible closures and maximum ties; independent horizontal slab x-interval merge and direct set difference. No producer imports, segment tree, radius/window truncation, or general C2 claim.')
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(dict(status=result['status'],sources=len(positions),maxima=len(winners),full_halo=full_area['display_only'],deltas=[r['actual_Delta']['display_only'] for r in jobs],seconds=result['elapsed_seconds'])))


if __name__=='__main__':main()
