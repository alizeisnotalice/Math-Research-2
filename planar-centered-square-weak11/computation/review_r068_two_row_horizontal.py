#!/usr/bin/env python3
"""Independent generic bbox-bitmaps + horizontal interval-union N128 replay.

No producer imports, segment tree, source-subset array, or rectangle-piece
subtraction. Pinned inputs use archive-relative paths; outputs use dataTemp.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse,gzip,hashlib,json,math,os,sys,time

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def need(v,msg):
    if not v:raise RuntimeError(msg)
def sqrtfloor(n):
    lo=0;hi=1<<((n.bit_length()+1)//2)
    while lo+1<hi:
        k=(lo+hi)//2
        if k*k<=n:lo=k
        else:hi=k
    return lo

def merge(seq):
    seq=sorted((a,b)for a,b in seq);out=[]
    for a,b in seq:
        if not out or a>out[-1][1]:out.append([a,b])
        elif b>out[-1][1]:out[-1][1]=b
    return out

def minus_interval(a,b,core):
    if core is None:return [(a,b)]
    c,d=core
    if d<=a or c>=b:return [(a,b)]
    return [(x,y)for x,y in ((a,min(b,c)),(max(a,d),b))if x<y]
def length(seq):return sum(b-a for a,b in seq)
def overlap(seq,core):
    if core is None:return 0
    c,d=core;return sum(max(0,min(b,d)-max(a,c))for a,b in seq)
def horizontal(internal,cross,Q):
    # cross records are (full rectangle, exact same-full-bbox residual core).
    ys=sorted({v for r in internal+[r for r,k in cross]+[k for r,k in cross if k is not None]+[Q]for v in r[2:]})
    sums={k:0 for k in ('H_S','H_full','W','Wprime','true_deletion','W_minus_Delta_overpayment','Wprime_minus_Delta_overpayment')}
    maxactive=0
    for y,z in zip(ys,ys[1:]):
        old=merge([(r[0],r[1])for r in internal if r[2]<=y<r[3]])
        cr=[];wr=[];qr=(Q[0],Q[1])if Q[2]<=y<Q[3]else None
        for r,k in cross:
            if not r[2]<=y<r[3]:continue
            cr.append((r[0],r[1]));kc=(k[0],k[1])if k is not None and k[2]<=y<k[3]else None
            wr.extend(minus_interval(r[0],r[1],kc))
        cr=merge(cr);wr=merge(old+wr)
        wp=merge(old+[piece for a,b in cr for piece in minus_interval(a,b,qr)])
        full=merge(old+cr+([]if qr is None else[qr]));width=z-y
        assert overlap(old,qr)==0
        wlen,plen=length(wr),length(wp);dlen=length(full)-(0 if qr is None else qr[1]-qr[0]);over=overlap(wr,qr);pover=overlap(wp,qr)
        assert plen==dlen and wlen==dlen+over and pover==0
        vals=dict(H_S=length(old),H_full=length(full),W=wlen,Wprime=plen,true_deletion=dlen,W_minus_Delta_overpayment=over,Wprime_minus_Delta_overpayment=pover)
        for k,v in vals.items():sums[k]+=v*width
        maxactive=max(maxactive,len(old)+len(cr))
    return sums,dict(y_strips=len(ys)-1,maximum_merged_old_plus_cross_intervals=maxactive)

def replay(data,cat,expected):
    N=data['N'];n=len(data['points']);U=32*N;assert N==128 and n==257
    p=[tuple(F(v)*U for v in q)for q in data['points']];assert all(v.denominator==1 for q in p for v in q);p=[tuple(map(int,q))for q in p]
    m=list(map(F,data['initial_masses']));assert m==[F(1,2*N)]*(2*N)+[F(1,4)];S=data['S_mask'];assert S==(1<<(2*N))-1
    # Generic coordinate-boundary intersections; no row/suffix descriptors.
    xs=sorted({x for x,y in p});ys=sorted({y for x,y in p});xgroups={x:sum(1<<i for i,(xx,yy)in enumerate(p)if xx==x)for x in xs};ygroups={y:sum(1<<i for i,(xx,yy)in enumerate(p)if yy==y)for y in ys}
    yranges=[]
    for i,c in enumerate(ys):
        mask=0
        for d in ys[i:]:mask|=ygroups[d];yranges.append((c,d,mask))
    boxes={}
    for i,a in enumerate(xs):
        xm=0
        for b in xs[i:]:
            xm|=xgroups[b]
            for c,d,ym in yranges:
                mask=xm&ym
                if not mask or any(not(mask&g)for g in (xgroups[a],xgroups[b],ygroups[c],ygroups[d])):continue
                key=(a,b,c,d);assert key not in boxes
                units=2*(mask&S).bit_count()+N*bool(mask&~S);span=max(b-a,d-c)
                boxes[key]=(mask,units,span*span<256*N*units)
    assert len(boxes)==expected['all_closed_bboxes']==25153
    assert cat['N']==N and cat['coordinate_unit']==f'1/{U}'and cat['mass_unit']==f'1/{4*N}'
    catalog={tuple(r[4:8]):r for r in cat['rows']};assert len(catalog)==len(cat['rows'])==len(boxes)and set(catalog)==set(boxes)
    D=int(expected['coordinate_scale_D']);assert D%U==0;mesh=D//U;assert mesh&(mesh-1)==0
    den=D*D;radii={};eligible=0;maximum=0;winners=[];classes={}
    internal=[[],[]];cross=[[],[]];qr=D//4;Q=(39*N*mesh-qr,39*N*mesh+qr,-qr,qr)
    for (a,b,c,d),(mask,units,ok)in boxes.items():
        row=catalog[(a,b,c,d)];assert row[3]==units and row[8]is ok
        if not(mask&S):kind,ii,jj=6,-1,-1
        elif mask&~S:
            kind={( -14*N,0):3,(0,14*N):4,(-14*N,14*N):5}[(c,d)];assert a%28==0;ii,jj=a//28,N-1
        else:
            kind={(-14*N,-14*N):0,(14*N,14*N):1,(-14*N,14*N):2}[(c,d)];assert a%28==b%28==0;ii,jj=a//28,b//28
        assert row[:3]==[kind,ii,jj]

        if not ok:assert len(row)==9;continue
        eligible+=1;kind=row[0];classes[str(kind)]=classes.get(str(kind),0)+1
        if units>maximum:maximum=units;winners=[mask]
        elif units==maximum:winners.append(mask)
        if not(mask&S):assert mask==1<<(2*N)and row[9]==qr;continue
        if units not in radii:radii[units]=sqrtfloor(64*N*units*mesh*mesh)
        k=radii[units];assert row[9]==k and k*k<=64*N*units*mesh*mesh<(k+1)**2
        K=None
        if mask&~S:
            core=(b*mesh-qr,a*mesh+qr,d*mesh-qr,c*mesh+qr)
            if core[0]<core[1]and core[2]<core[3]:K=core
        for side,r in enumerate((max(0,k-1),k+1)):
            R=(b*mesh-r,a*mesh+r,d*mesh-r,c*mesh+r)
            if R[0]>=R[1]or R[2]>=R[3]:continue
            if mask&~S:cross[side].append((R,K))
            else:internal[side].append(R)
    assert eligible==expected['strict_eligible_bboxes']==15205 and maximum==4*N and winners==expected['maximizer_masks']==[S]
    assert classes==expected['eligible_by_class']and len(radii)==expected['distinct_mass_radii']==139
    values=[];stats=[]
    for side in (0,1):
        a,s=horizontal(internal[side],cross[side],Q);values.append(a);stats.append(s)
    metrics={k:dict(lower=str(F(values[0][k],den)),upper=str(F(values[1][k],den)))for k in values[0]}
    for k,v in metrics.items():assert expected[k]==v,(k,v,expected[k])
    assert expected['remaining_halo_area']=='1/4'and expected['Wprime_equals_Delta_at_both_directed_levels']is True
    assert max(r[1]for r in internal[1])<Q[0]
    # Family carrier numerical sign only; the analytic covering argument is
    # reviewed separately, not inferred from this one N computation.
    H=1<<128;k2=sqrtfloor(2*H*H);k3=sqrtfloor(3*H*H)
    carrier_hi=F(399,128)+F(7,8)*F(k2+1,H)-F(25,32)*F(k3,H)
    assert carrier_hi<F(389,128)<F(7,2)
    return dict(N=N,source_count=n,generic_closed_bboxes=len(boxes),eligible_bboxes=eligible,all_maximum_masks=winners,distinct_mass_radii=len(radii),metrics=metrics,horizontal_statistics=stats,general_family_carrier_sign_upper=str(carrier_hi),all_metric_endpoints_exact_match=True,scope='One actual N128 independent horizontal-stripe replay and generic coordinate-boundary catalog census. No segmenttree/producer imports or source-subset enumeration. All-N carrier geometry is a separate analytic proof obligation.')
def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--archive-root',type=Path,default=Path(__file__).resolve().parents[2]);ap.add_argument('--spec',type=Path,default=Path(__file__).with_name('r068-two-row-replay-inputs.json'));ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
    need(sys.flags.optimize==0 and Path('/Volumes/data').is_mount(),'Assertions/data volume required');root=args.archive_root.resolve()
    if root.name=='planar-centered-square-weak11':root=root.parent
    out=args.output.resolve();rp=out.with_suffix('.receipt.json');need(out.is_relative_to(Path('/Volumes/data/CodexTemp').resolve())and out.suffix=='.json'and not out.exists()and not rp.exists(),'New dataTemp JSON required');out.parent.mkdir(parents=True,exist_ok=True)
    start=time.perf_counter();rec=dict(pid=os.getpid(),python=sys.version,entry_sha256=sha(Path(__file__)),spec_sha256=sha(args.spec),archive_root=str(root),command=sys.argv)
    try:
        spec=json.loads(args.spec.read_text())
        for rel,h in spec['dependencies'].items():
            p=(root/rel).resolve();need(p.is_relative_to(root)and p.is_file()and sha(p)==h,f'Archive dependency mismatch {rel}')
        with gzip.open(root/spec['catalog'],'rt')as f:catalog=json.load(f)
        expected=json.loads((root/spec['expected_result']).read_text())
        need(expected['source_sha256']=='b0c9bc02ea6219bb6c20a6a8bb0b3273fc0fb942859ea2f46318a3a2952e3efc','Frozen producer source pin')
        need(expected['input_sha256']==sha(root/spec['input'])and expected['catalog_sha256']==sha(root/spec['catalog']),'N128 input/catalog output chain')
        result=replay(json.loads((root/spec['input']).read_text()),catalog,expected)
        result.update(status='PASS_R68_N128_INDEPENDENT_HORIZONTAL_UNION',**rec,elapsed_seconds=time.perf_counter()-start,dependency_pins=spec['dependencies']);out.write_text(json.dumps(result,indent=2)+'\n')
        rec.update(status=result['status'],elapsed_seconds=result['elapsed_seconds'],output_sha256=sha(out));rp.write_text(json.dumps(rec,indent=2)+'\n');print(json.dumps(rec,indent=2))
    except BaseException as e:
        rec.update(status='FAIL_R68_N128_REPLAY',error=str(e));rp.write_text(json.dumps(rec,indent=2)+'\n');raise
if __name__=='__main__':main()
