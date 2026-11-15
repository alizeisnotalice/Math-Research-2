#!/usr/bin/env python3
"""Exact N=128 repeated-capacity family, O(N^2 log N) union sweep.

All closed source bboxes are contiguous old intervals in one/both rows,
or the external point plus a suffix in one/both rows, or the singleton.
No source-subset enumeration or general 2^n cut arrays are used.
"""
import argparse
from collections import Counter
from fractions import Fraction as F
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import sys
import time


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def union_area(rects):
    """Exact vertical segment-tree coverage swept across integer x events."""
    if not rects:return 0,{'rectangles':0,'y_intervals':0,'x_events':0}
    ys=sorted({r[j] for r in rects for j in (2,3)});ranks={v:i for i,v in enumerate(ys)}
    events=[]
    for x0,x1,y0,y1 in rects:
        assert x0<x1 and y0<y1
        events.append((x0,1,ranks[y0],ranks[y1]));events.append((x1,-1,ranks[y0],ranks[y1]))
    events.sort();M=len(ys)-1;cover=[0]*(4*M+8);length=[0]*(4*M+8)
    def update(node,left,right,a,b,delta):
        if a<=left and right<=b:
            cover[node]+=delta;assert cover[node]>=0
        else:
            middle=(left+right)//2
            if a<middle:update(node*2,left,middle,a,b,delta)
            if middle<b:update(node*2+1,middle,right,a,b,delta)
        if cover[node]:length[node]=ys[right]-ys[left]
        elif right-left==1:length[node]=0
        else:length[node]=length[node*2]+length[node*2+1]
    area=0;previous=events[0][0];j=0;xs=0
    while j<len(events):
        x=events[j][0];area+=(x-previous)*length[1];xs+=1
        while j<len(events) and events[j][0]==x:
            _,delta,a,b=events[j];update(1,0,M,a,b,delta);j+=1
        previous=x
    assert length[1]==0 and not any(cover)
    return area,{'rectangles':len(rects),'y_intervals':M,'x_events':xs}


def clip(rect,Q):
    r=(max(rect[0],Q[0]),min(rect[1],Q[1]),max(rect[2],Q[2]),min(rect[3],Q[3]))
    return r if r[0]<r[1] and r[2]<r[3] else None


def difference(rect,Q):
    if Q is None:return [rect]
    I=clip(rect,Q)
    if I is None:return [rect]
    a,b,c,d=rect;x0,x1,y0,y1=I
    pieces=[(a,x0,c,d),(x1,b,c,d),(x0,x1,c,y0),(x0,x1,y1,d)]
    return [r for r in pieces if r[0]<r[1] and r[2]<r[3]]


def enclose(lo,hi,den):
    assert 0<=lo<=hi
    return {'lower':str(F(lo,den)),'upper':str(F(hi,den))}


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--N',type=int,default=128)
    ap.add_argument('--bits',type=int,default=64);ap.add_argument('--out-dir',type=Path,required=True)
    args=ap.parse_args();assert os.path.ismount('/Volumes/data')
    N=args.N;bits=args.bits;assert N>=32 and N%32==0 and bits>=32
    started=time.perf_counter();U=32*N;Qmesh=1<<bits;D=U*Qmesh;den=D*D
    # Source coordinates are integers / U. Each old mass has units2 and the
    # external mass has unitsN, where one mass unit is 1/(4N).
    points=[(28*i,-14*N) for i in range(N)]+[(28*i,14*N) for i in range(N)]+[(39*N,0)]
    masses=[F(1,2*N)]*(2*N)+[F(1,4)];S=(1<<(2*N))-1;external=1<<(2*N)
    source={'schema':'r068-two-row-repeated-capacity-input-v1','N':N,'source_count':2*N+1,
            'points':[[str(F(x,U)),str(F(y,U))] for x,y in points],
            'initial_masses':list(map(str,masses)),'S_mask':S,'S_indices':list(range(2*N)),
            'family':'p_i^+-=(7i/(8N),+-7/16),mass1/(2N);q=(39/32,0),mass1/4',
            'scope':'One prescribed N, no random inputs; source catalog is enumerated by continuous intervals and suffixes.'}
    # kind0/1 lower/upper old interval, kind2 both; kind3/4 lower/upper
    # external suffix, kind5 both; kind6 external singleton.
    descriptors=[]
    for i in range(N):
        for j in range(i,N):
            ell=j-i+1;a=28*i;b=28*j;row=((1<<ell)-1)<<i
            descriptors.extend([(0,i,j,2*ell,a,b,-14*N,-14*N,row),
                                (1,i,j,2*ell,a,b,14*N,14*N,row<<N),
                                (2,i,j,4*ell,a,b,-14*N,14*N,row|(row<<N))])
    for k in range(N):
        ell=N-k;a=28*k;b=39*N;row=((1<<ell)-1)<<k
        descriptors.extend([(3,k,N-1,N+2*ell,a,b,-14*N,0,row|external),
                            (4,k,N-1,N+2*ell,a,b,0,14*N,(row<<N)|external),
                            (5,k,N-1,N+4*ell,a,b,-14*N,14*N,row|(row<<N)|external)])
    descriptors.append((6,-1,-1,N,39*N,39*N,0,0,external))
    assert len(descriptors)==3*N*(N+1)//2+3*N+1
    assert len({r[-1] for r in descriptors})==len(descriptors)
    radcache={}
    def radical(units):
        if units not in radcache:
            square=64*N*units*Qmesh*Qmesh;k=math.isqrt(square)
            assert k*k<=square<(k+1)*(k+1)
            radcache[units]=(k,max(0,k-1),k+1)
        return radcache[units]
    internal=[[],[]];cross=[[],[]];W=[[],[]];WP=[[],[]]
    eligible=[];best=0;maxima=[];types=Counter();catalog=[]
    qradius=8*N*Qmesh;Q=(39*N*Qmesh-qradius,39*N*Qmesh+qradius,-qradius,qradius)
    assert (Q[1]-Q[0])*(Q[3]-Q[2])==den//4
    empty_old_kernel=0;nonempty_old_kernel=0
    for kind,i,j,units,a,b,c,d,mask in descriptors:
        span=max(b-a,d-c);ok=span*span<256*N*units
        row=[kind,i,j,units,a,b,c,d,ok]
        if not ok:catalog.append(row);continue
        types[kind]+=1;eligible.append((kind,i,j,units,mask))
        if units>best:best=units;maxima=[mask]
        elif units==best:maxima.append(mask)
        if kind==6:
            row.append(qradius);catalog.append(row);continue
        k,ri,ro=radical(units);row.append(k);catalog.append(row)
        kernel=None
        if kind>=3:
            # b=1/4. Same-full-bbox residual core has exact rational radius.
            K=(b*Qmesh-qradius,a*Qmesh+qradius,d*Qmesh-qradius,c*Qmesh+qradius)
            if K[0]<K[1] and K[2]<K[3]:
                kernel=K;nonempty_old_kernel+=1
                assert Q[0]<=K[0]<=K[1]<=Q[1] and Q[2]<=K[2]<=K[3]<=Q[3]
            else:empty_old_kernel+=1
        for side,radius in enumerate((ri,ro)):
            R=(b*Qmesh-radius,a*Qmesh+radius,d*Qmesh-radius,c*Qmesh+radius)
            if not(R[0]<R[1] and R[2]<R[3]):continue
            if kind<=2:
                internal[side].append(R);W[side].append(R);WP[side].append(R)
            else:
                cross[side].append(R);W[side]+=difference(R,kernel);WP[side]+=difference(R,Q)
    assert best==4*N and maxima==[S]
    assert 16*(28*(N-1))**2>9*U*U and 16*(28*N)**2>9*U*U
    assert F(39,32)**2>F(35,33)*F(5,4)
    selected=range(N//4+1,15*N//32)
    selected_mass_sum=sum((F(N-k,N) for k in selected),F(0))
    assert len(selected)==7*N//32-1 and selected_mass_sum==F(41,64)*(7*N//32-1)
    # Exact enclosure-level separation, consistent with the general proof
    # H(S).right<=53/56<31/32=Q.left.
    assert max(r[1] for r in internal[1])<Q[0]
    full=[internal[s]+cross[s]+[Q] for s in (0,1)]
    areas={};sweeps={}
    for name,lists in [('old_halo',internal),('full_halo',full),('W',W),('Wprime',WP)]:
        values=[]
        for side,rects in enumerate(lists):
            t=time.perf_counter();value,stats=union_area(rects);stats['seconds']=time.perf_counter()-t
            values.append(value);sweeps[f'{name}_{side}']=stats
        areas[name]=values
    residual_area=den//4
    delta=[v-residual_area for v in areas['full_halo']]
    assert areas['Wprime']==delta
    assert all(areas['Wprime'][s]<=areas['W'][s]<=areas['full_halo'][s] for s in (0,1))
    over=[]
    for side in (0,1):
        rects=[c for r in W[side] if (c:=clip(r,Q)) is not None]
        value,stats=union_area(rects);over.append(value);sweeps[f'W_overlap_Q_{side}']=stats
    assert all(areas['W'][s]==delta[s]+over[s] for s in (0,1))
    # A global six-carrier bound for this family, independent of N.
    H=1<<128;k2=math.isqrt(2*H*H);k3=math.isqrt(3*H*H)
    r2lo,r2hi=F(k2,H),F(k2+1,H);r3lo,r3hi=F(k3,H),F(k3+1,H)
    carrier_lo=F(399,128)+F(7,8)*r2lo-F(25,32)*r3hi
    carrier_hi=F(399,128)+F(7,8)*r2hi-F(25,32)*r3lo
    assert carrier_hi<F(389,128)<F(7,2)
    args.out_dir.mkdir(parents=True,exist_ok=True)
    ip=args.out_dir/'input.json';ip.write_text(json.dumps(source,indent=2)+'\n')
    cp=args.out_dir/'catalog.json.gz'
    with cp.open('wb') as raw:
        with gzip.GzipFile(filename='',mode='wb',fileobj=raw,mtime=0) as z:
            z.write(json.dumps({'schema':'r068-complete-special-bbox-catalog-v1',
                'row_fields':['kind','i','j','mass_units','a','b','c','d','strict_eligible','radius_floor_if_eligible'],
                'coordinate_unit':f'1/{U}','mass_unit':f'1/{4*N}','N':N,'rows':catalog},separators=(',',':')).encode())
    payload={'schema':'r068-two-row-exact-erosion-union-v1','status':'PASS','N':N,'source_count':2*N+1,
             'python':sys.version,'platform':platform.platform(),'source_sha256':sha(Path(__file__)),
             'input_sha256':sha(ip),'catalog_sha256':sha(cp),'all_closed_bboxes':len(catalog),
             'strict_eligible_bboxes':len(eligible),'eligible_by_class':dict(types),
             'maximum_eligible_mass':'1','maximizer_masks':maxima,'S_rank':2*N,
             'S_widths':[str(F(7*(N-1),8*N)),'7/8'],'whole_mass':'5/4',
             'connected':True,'connected_witness':'eligible fullS plus any selected external both-row suffix',
             'nonterminal_margin':str(F(39,32)**2-F(35,33)*F(5,4)),
             'coordinate_scale_D':str(D),'root_grid_step':str(F(1,D)),
             'all_radius_floor_squared_checks_passed':True,'distinct_mass_radii':len(radcache),
             'W':enclose(*areas['W'],den),'Wprime':enclose(*areas['Wprime'],den),
             'true_deletion':enclose(*delta,den),'H_S':enclose(*areas['old_halo'],den),
             'H_full':enclose(*areas['full_halo'],den),'remaining_halo_area':'1/4',
             'W_minus_Delta_overpayment':enclose(*over,den),'Wprime_minus_Delta_overpayment':{'lower':'0','upper':'0'},
             'Wprime_equals_Delta_at_both_directed_levels':True,
             'old_halo_separation':{'analytic_old_right_upper':'53/56','Q_left':'31/32','strict_gap':'5/224',
                                    'directed_outer_old_right':str(F(max(r[1] for r in internal[1]),D))},
             'selected_repeated_capacities':{'count':len(selected),'sum_removed_masses':str(selected_mass_sum)},
             'crossing_full_bbox_kernel_counts':{'empty':empty_old_kernel,'nonempty':nonempty_old_kernel},
             'W_certifies_7over2':F(areas['W'][1],den)<F(7,2),
             'Wprime_certifies_7over2':F(areas['Wprime'][1],den)<F(7,2),
             'general_family_six_carrier_bound':{'expression':'399/128+7sqrt(2)/8-25sqrt(3)/32',
                 'lower':str(carrier_lo),'upper':str(carrier_hi),'rational_strict_upper':'389/128'},
             'sweep_statistics':sweeps,'wall_seconds':time.perf_counter()-started,
             'scope':f'One actual N={N}-family area evaluation. Complete continuous-interval/tail catalog; segment-tree union retains all overlaps. No2^source enumeration, no new general C2 bound.'}
    (args.out_dir/'result.json').write_text(json.dumps(payload,indent=2)+'\n')
    print(json.dumps({'status':'PASS','N':N,'closed':len(catalog),'eligible':len(eligible),
                      'W_upper_display':float(F(areas['W'][1],den)),
                      'Wprime_upper_display':float(F(areas['Wprime'][1],den)),
                      'overpayment_lower_display':float(F(over[0],den)),'seconds':payload['wall_seconds']}))


if __name__=='__main__':main()
