#!/usr/bin/env python3
"""Portable independent 2D-grid replay of pinned R67 W/Wprime states.

Standard library only. All read dependencies are pinned archive-relative paths.
--output names a NEW JSON file below /Volumes/data/CodexTemp.
"""
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction as F
from pathlib import Path
import argparse,hashlib,importlib.util,json,os,sys,time
sys.dont_write_bytecode=True
CORE='planar-centered-square-weak11/computation/review_r066_grid_budgets.py'
HELPER='planar-centered-square-weak11/computation/review_r064_mass_variation.py'
CORE_SHA='c6524c436d8efd15a1b941c3e37448608622f6723fc51e3325a6eaa44aa159b1'
HELPER_SHA='267b430e167801314bba81200250945bedcb1473918992f57ac03c1feb94cad4'
STATE_SHA='ab0e7dba2993869d69871384e9f95f0f095563b453dbde897561097d06251ccf'
GEOM_SHA='4ef61b815b9c4a9a8972cdaf31716ae7917b59243a29fb8697c8042137788246'
CLIP_SHA='4dc4fb494736c40334359f9950f052ee2e4a53b503883ffbc3f9bc0da69b1361'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def require(v,msg):
    if not v:raise RuntimeError(msg)
def boot(path):
    global A
    s=importlib.util.spec_from_file_location('own_frozen_grid',path);A=importlib.util.module_from_spec(s);s.loader.exec_module(A);A.init_worker()
def grid(layers,nlabels):
    xs=sorted({r[k]for lay in layers for r,w in lay for k in (0,1)})
    ys=sorted({r[k]for lay in layers for r,w in lay for k in (2,3)})
    nx,ny=len(xs),len(ys);require(nx*ny<4000000,'Grid memory guard')
    xi={v:i for i,v in enumerate(xs)};yi={v:i for i,v in enumerate(ys)}
    ar=[[0]*(nx*ny)for _ in layers]
    for a,lay in zip(ar,layers):
        for r,w in lay:
            x,z,y,t=xi[r[0]],xi[r[1]],yi[r[2]],yi[r[3]]
            for xx,yy,c in ((x,y,w),(z,y,-w),(x,t,-w),(z,t,w)):a[xx*ny+yy]+=c
    totals=[0]*16;ind=[[0]*nlabels for _ in range(4)];zero=[0]*6
    for x in range(nx):
        for y in range(ny):
            j=x*ny+y
            for a in ar:
                if x:a[j]+=a[j-ny]
                if y:a[j]+=a[j-1]
                if x and y:a[j]-=a[j-ny-1]
                assert a[j]>=0
            if x+1==nx or y+1==ny:continue
            fi,fo,ri,ro,mi,mo,bi,bo,ti,to,zi,zo=[a[j]for a in ar]
            wi,wo,pi,po=mi&~bo,mo&~bi,mi&~to,mo&~ti
            area=(xs[x+1]-xs[x])*(ys[y+1]-ys[y])
            flags=(wi,wo,pi,po,fi and not ro,fo and not ri,fi,fo,ri,ro,wi and ri,wo and ro,pi and ri,po and ro,wi and not po,wo and not pi)
            for k,v in enumerate(flags):
                if v:totals[k]+=area
            checks=(pi and not wi,po and not wo,wo and not fo,fi and not ro and not pi,wo and not zo,pi and not zo)
            for k,v in enumerate(checks):
                if v:zero[k]+=area
            for dest,mask in zip(ind,(wi,wo,pi,po)):
                while mask:
                    bit=mask&-mask;dest[bit.bit_length()-1]+=area;mask^=bit
    assert not any(zero),zero
    return totals,ind,zero,(nx-1)*(ny-1)
def evaluate(job):
    data,saved=job;n=len(data['points']);p=[tuple(map(F,v))for v in data['points']];m=list(map(F,data['initial_masses']))
    assert len(set(p))==n and len(m)==n and min(m)>0
    S=data['S_mask'];assert 0<S<(1<<n)-1
    boxes=A.G.boxes(p);mass={s:sum(m[i]for i in b['ids'])for s,b in boxes.items()}
    elig=[s for s,b in boxes.items()if mass[s]>b['d2']];Lmax=max(mass[s]for s in elig)
    winners=sorted(s for s in elig if mass[s]==Lmax);L=sum(m[i]for i in range(n)if S>>i&1)
    assert saved['source_count']==n and saved['S_mask']==S and saved['S_rank']==S.bit_count()
    assert F(saved['L'])==F(saved['maximum_eligible_mass'])==L
    assert L==Lmax and S in winners and saved['all_maximum_masks']==winners
    assert saved['closed_bbox_count']==len(boxes)and saved['strict_eligible_count']==len(elig)
    bits=saved['bits'];D2=1<<(2*bits);labels=sorted(s for s in elig if s&S)
    assert saved['intersect_S_eligible_count']==len(labels)
    full={s:A.rect(boxes[s]['bbox'],mass[s],bits)for s in elig}
    keep=[i for i in range(n)if not S>>i&1];rp=[p[i]for i in keep];rm=[m[i]for i in keep];rb=A.G.boxes(rp)
    remain=[A.rect(b['bbox'],sum(rm[i]for i in b['ids']),bits)for b in rb.values()]
    layers=[[]for _ in range(12)]
    for rr in full.values():
        for k in (0,1):
            if rr[k]:layers[k].append((rr[k],1))
    for rr in remain:
        for k in (0,1):
            if rr[k]:layers[2+k].append((rr[k],1))
    labelmeta=[]
    for j,s in enumerate(labels):
        b=boxes[s];ids=[i for i in b['ids']if not S>>i&1];bm=sum((m[i]for i in ids),F(0));am=mass[s]-bm
        core=A.rect(b['bbox'],bm,bits)
        if ids:
            tight=(min(p[i][0]for i in ids),max(p[i][0]for i in ids),min(p[i][1]for i in ids),max(p[i][1]for i in ids))
            tightcore=A.rect(tight,bm,bits);td2=max(tight[1]-tight[0],tight[3]-tight[2])**2
        else:tightcore=[None,None];td2=None
        for k in (0,1):
            for index,rr in ((4,full[s]),(6,core),(8,tightcore)):
                if rr[k]:layers[index+k].append((rr[k],1<<j))
        labelmeta.append(dict(mask=s,rank=s.bit_count(),intersect_S_rank=(s&S).bit_count(),mass=str(mass[s]),a_intersection_mass=str(am),b_outside_mass=str(bm),full_bbox_span_squared=str(b['d2']),outside_tight_bbox_span_squared=None if td2 is None else str(td2),full_bbox_kernel_nonempty=b['d2']<bm,tight_outside_kernel_nonempty=bool(ids)and td2<bm,analytic_single_label_upper_a=str(am)))
    for i in range(n):
        if S>>i&1:
            rr=A.rect((p[i][0],p[i][0],p[i][1],p[i][1]),L,bits)
            for k in (0,1):
                if rr[k]:layers[10+k].append((rr[k],1))
    vals,ind,zero,cells=grid(layers,len(labels))
    def interval(lo,hi):
        assert 0<=lo<=hi
        return dict(lower=str(F(lo,D2)),upper=str(F(hi,D2)),normalized_lower=str(F(lo,D2)/L),normalized_upper=str(F(hi,D2)/L))
    names=('W','Wprime','true_deletion','full_halo','remaining_halo','W_minus_Delta_overpayment','Wprime_minus_Delta_overpayment','W_minus_Wprime_refinement')
    metrics={name:interval(vals[2*j],vals[2*j+1])for j,name in enumerate(names)}
    for k,v in metrics.items():assert saved[k]==v,(data['id'],k)
    slabels={r['mask']:r for r in saved['all_label_rings']};assert set(slabels)==set(labels)
    for j,row in enumerate(labelmeta):
        row['W_label_ring']=interval(ind[0][j],ind[1][j]);row['Wprime_label_ring']=interval(ind[2][j],ind[3][j])
        assert row==slabels[row['mask']],(data['id'],row['mask'])
        assert F(ind[0][j],D2)<=F(row['a_intersection_mass'])and F(ind[2][j],D2)<=F(row['a_intersection_mass'])
    for w,offset in (('W',0),('Wprime',2)):
        residual=[str(F(vals[offset]-vals[5]-vals[11 if w=='W'else 13],D2)),str(F(vals[offset+1]-vals[4]-vals[10 if w=='W'else 12],D2))]
        assert list(map(F,residual))[0]<=0<=list(map(F,residual))[1]
        assert saved['overpayment_identity_'+w+'_residual']==residual
        assert saved['certifies_'+w+'_le_7over2_L']==(F(vals[offset+1],D2)<=F(7,2)*L)
    assert all(F(v)==0 for v in saved['exact_zero_area_inclusion_checks'].values())
    comps=A.components(p,m,elig,boxes);connected=len(comps)==1
    ss=boxes[S]['bbox'];fat=16*(ss[1]-ss[0])**2>9*L and 16*(ss[3]-ss[2])**2>9*L
    padding=None
    if data['id']=='r067-positive-nine-plus-one-uniform-neighborhood-obstruction':
        d,eta,eps=F(1023,1024),F(1,1048576),F(1,1024)
        gridpoints={(x,y)for x in (F(0),d/2,d)for y in (F(0),d/2,d)}
        assert {p[i]for i in range(n)if S>>i&1}==gridpoints and S==511 and n==10
        assert sum(m[i]for i in range(9))==1 and m[0]==F(1,5)and m[8]==F(4,5)-7*eta and m[9]==eps
        subsets=[]
        for mask in range(1,1<<n):
            ids=[i for i in range(n)if mask>>i&1];mm=sum(m[i]for i in ids)
            if mm>A.G.span2(p,ids):subsets.append((mask,mm))
        assert max(mm for mask,mm in subsets)==1 and [mask for mask,mm in subsets if mm==1]==[511]
        assert (d/2)**2<F(3,4)
        excess=F(246397,131072000)-21*eta/4
        assert excess==F(982963,524288000)>0
        # All full-bbox residual cores vanish; residual singleton is covered
        # by heavy internal label, so W=H and Wprime=W for this input.
        assert all(not r['full_bbox_kernel_nonempty']for r in labelmeta)
        assert metrics['W']==metrics['Wprime']==metrics['full_halo']
        padding=dict(nonempty_subsets=1023,support_cartesian_overlap=True,filled_support_budget_excess_lower=str(excess),analytic_W_equals_H=True)
    return dict(id=data['id'],sources=n,S_mask=S,S_rank=S.bit_count(),L=str(L),maximum_masks=winners,connected=connected,whole_nonterminal=33*A.G.span2(p,range(n))>35*sum(m),both_axes_fat=fat,labels=len(labels),metrics=metrics,zero_inclusion_areas=zero,compressed_grid_cells=cells,all_label_rows=labelmeta,positive_padding_checks=padding)
def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--archive-root',type=Path,default=Path(__file__).resolve().parents[2]);ap.add_argument('--spec',type=Path,default=Path(__file__).with_name('r067-independent-replay-inputs.json'));ap.add_argument('--output',type=Path,required=True);ap.add_argument('--workers',type=int,default=2)
    args=ap.parse_args();require(sys.flags.optimize==0,'Python -O unsupported');require(Path('/Volumes/data').is_mount(),'data not mounted');require(1<=args.workers<=12,'workers1..12')
    root=args.archive_root.resolve()
    if root.name=='planar-centered-square-weak11':root=root.parent
    out=args.output.expanduser().resolve();require(out.is_relative_to(Path('/Volumes/data/CodexTemp').resolve())and out.suffix=='.json','New output must be a JSON file below dataTemp')
    receipt=out.with_suffix('.receipt.json');require(not out.exists()and not receipt.exists(),'Choose new output paths');out.parent.mkdir(parents=True,exist_ok=True)
    start=time.perf_counter();meta=dict(pid=os.getpid(),python=sys.version,entry_sha256=sha(Path(__file__)),spec_sha256=sha(args.spec),archive_root=str(root),command=sys.argv,workers=args.workers)
    try:
        spec=json.loads(args.spec.read_text());pins=spec['dependencies'];require(pins[CORE]==CORE_SHA and pins[HELPER]==HELPER_SHA,'Frozen core dependency pins')
        for rel,h in pins.items():
            q=(root/rel).resolve();require(q.is_relative_to(root)and q.is_file(),f'Missing archive dependency {rel}');require(sha(q)==h,f'Dependency SHA mismatch {rel}')
        jobs=[]
        for c in spec['cases']:
            require(c['input']in pins and c['result']in pins,'Unpinned input/result')
            ip,op=root/c['input'],root/c['result'];d=json.loads(ip.read_text());r=json.loads(op.read_text())
            require(r['source_sha256']==STATE_SHA and r['geometry_source_sha256']==GEOM_SHA and r['rectclip_source_sha256']==CLIP_SHA,'Production source chain')
            require(r['input_sha256']==sha(ip)and r['status']=='PASS','Input/output chain')
            jobs.append((d,r))
        core=str(root/CORE)
        if args.workers==1:boot(core);rows=[evaluate(j)for j in jobs]
        else:
            with ProcessPoolExecutor(max_workers=args.workers,initializer=boot,initargs=(core,))as pool:rows=list(pool.map(evaluate,jobs))
        require(F(7095,2048)==(1+F(217,256))*(1+F(7,8))<F(7,2),'N32 analytic container bound')
        result=dict(status='PASS_INDEPENDENT_R67_2D_EROSION_STATES',**meta,elapsed_seconds=time.perf_counter()-start,cases=rows,statistics=dict(cases=len(rows),labels=sum(r['labels']for r in rows),compressed_grid_cells=sum(r['compressed_grid_cells']for r in rows)),dependency_pins=pins,scope='All specified labels, actual full/remaining halos, W/Wprime and old-halo overpayment rebuilt with independent coordinate2D difference/bitset algorithm. No producer strip/rect-difference import. Finite input certificates only.')
        out.write_text(json.dumps(result,indent=2)+'\n');meta.update(status=result['status'],elapsed_seconds=result['elapsed_seconds'],output_sha256=sha(out),statistics=result['statistics']);receipt.write_text(json.dumps(meta,indent=2)+'\n');print(json.dumps(meta,indent=2))
    except BaseException as e:
        meta.update(status='FAIL_R67_INDEPENDENT_REPLAY',error=str(e),elapsed_seconds=time.perf_counter()-start);receipt.write_text(json.dumps(meta,indent=2)+'\n');raise
if __name__=='__main__':main()
