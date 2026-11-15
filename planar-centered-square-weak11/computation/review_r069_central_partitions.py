#!/usr/bin/env python3
"""Independent14-layer2D source/area/K-partition replay, no producer import.

Archive-relative pinned inputs; mandatory new JSON output below dataTemp.
"""
from pathlib import Path
from fractions import Fraction as F
from concurrent.futures import ProcessPoolExecutor
import argparse,hashlib,importlib.util,json,math,os,sys,time
sys.dont_write_bytecode=True
PREFIX='planar-centered-square-weak11/computation/'
CORE_PINS={PREFIX+'review_r068_state_geometry.py':'889b8132e92f6c48cdf1537857ee0cd4a6440b94aa6571d8fc9135cfc252d822',PREFIX+'review_r067_erosion_states.py':'513e59dad17949c419ba98ec9c0e9ef3738871b92f32074df35eed874fd48eb0',PREFIX+'review_r066_grid_budgets.py':'c6524c436d8efd15a1b941c3e37448608622f6723fc51e3325a6eaa44aa159b1',PREFIX+'review_r064_mass_variation.py':'267b430e167801314bba81200250945bedcb1473918992f57ac03c1feb94cad4'}
PRODUCER_SHA='c7f64b45aa99f19e482e70479288c6e6e0bcaae0c780f35e983835630fba5ac7'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def need(v,msg):
    if not v:raise RuntimeError(msg)
def boot(root):
    global C,A
    s=importlib.util.spec_from_file_location('own_frozen_r068_geometry',Path(root)/(PREFIX+'review_r068_state_geometry.py'));C=importlib.util.module_from_spec(s);s.loader.exec_module(C);C.boot(Path(root));A=C.A

def area_grid(layers):
    xs=sorted({r[k]for lay in layers for r,w in lay for k in (0,1)});ys=sorted({r[k]for lay in layers for r,w in lay for k in (2,3)})
    nx,ny=len(xs),len(ys);need(nx*ny<2000000,'Grid memory guard');xi={v:i for i,v in enumerate(xs)};yi={v:i for i,v in enumerate(ys)}
    aa=[[0]*(nx*ny)for _ in layers]
    for ar,lay in zip(aa,layers):
        for r,w in lay:
            xl,xr,yl,yr=xi[r[0]],xi[r[1]],yi[r[2]],yi[r[3]]
            for x,y,v in ((xl,yl,w),(xr,yl,-w),(xl,yr,-w),(xr,yr,w)):ar[x*ny+y]+=v
    names=('Hfull','R','I','W','Wprime','Delta','P_pure_mixed','internal_old_credit','W_old_surplus','Wprime_old_surplus','Wprime_surplus_outside_I');totals={name:[[0,0]for _ in range(3)]for name in names};zero=[0]*6
    for x in range(nx):
        for y in range(ny):
            j=x*ny+y
            for ar in aa:
                if x:ar[j]+=ar[j-ny]
                if y:ar[j]+=ar[j-1]
                if x and y:ar[j]-=ar[j-ny-1]
                assert ar[j]>=0
            if x+1==nx or y+1==ny:continue
            fi,fo,ri,ro,ii,io,mi,mo,bi,bo,ti,to,ki,ko=[ar[j]for ar in aa]
            wi,wo,pi,po=mi&~bo,mo&~bi,mi&~to,mo&~ti
            flags=((fi,fo),(ri,ro),(ii,io),(wi,wo),(pi,po),(fi and not ro,fo and not ri),(fi and not io and not ro,fo and not ii and not ri),(ii and ri,io and ro),(wi and ri,wo and ro),(pi and ri,po and ro),(pi and ri and not io,po and ro and not ii))
            area=(xs[x+1]-xs[x])*(ys[y+1]-ys[y])
            for name,(lo,hi)in zip(names,flags):
                assert not lo or hi
                t=totals[name]
                if lo:t[0][0]+=area
                if hi:t[0][1]+=area
                if lo and ki:t[1][0]+=area
                if hi and ko:t[1][1]+=area
                if lo and not ko:t[2][0]+=area
                if hi and not ki:t[2][1]+=area
            checks=(io and not fo,ro and not fo,pi and not wi,po and not wo,wo and not fo,fi and not ro and not pi)
            for k,v in enumerate(checks):
                if v:zero[k]+=area
    assert not any(zero),zero
    return totals,zero,(nx-1)*(ny-1)

def evaluate(job):
    data,saved=job;start=time.perf_counter();p=[tuple(map(F,v))for v in data['points']];m=list(map(F,data['initial_masses']));n=len(p);S=data['S_mask']
    assert len(m)==n and min(m)>0 and len(set(p))==n and 0<S<(1<<n)-1
    assert saved['id']==data['id']and saved['source_input_id']==data['source_input_id']and saved['source_count']==n and saved['S_mask']==S and saved['S_rank']==S.bit_count()
    boxes=A.G.boxes(p);mm={s:sum(m[i]for i in b['ids'])for s,b in boxes.items()};elig=[s for s,b in boxes.items()if mm[s]>b['d2']]
    L=max(mm[s]for s in elig);winners=sorted(s for s in elig if mm[s]==L)
    assert S in winners and winners==saved['all_maximum_masks']==data['expected_maximum_masks']and sum(m[i]for i in range(n)if S>>i&1)==L==F(data['expected_L'])==F(saved['L'])
    c=math.lcm(*(v.denominator for q in p for v in q));bits=saved['bits_beyond_rational_source_grid'];D=c*(1<<bits);D2=D*D
    assert c==saved['computational_integer_dilation']and F(saved['physical_grid_step'])==F(1,D)
    labels=sorted(s for s in elig if s&S);full={s:C.receiver(boxes[s]['bbox'],mm[s],D)for s in elig};layers=[[]for _ in range(14)]
    for rr in full.values():
        for k in (0,1):
            if rr[k]:layers[k].append((rr[k],1))
    ids=[i for i in range(n)if S>>i&1];keep=[i for i in range(n)if not S>>i&1]
    ins=[p[i]for i in ids];im=[m[i]for i in ids];rem=[p[i]for i in keep];rm=[m[i]for i in keep];ib=A.G.boxes(ins);rb=A.G.boxes(rem)
    for offset,points,mass,bb in ((2,rem,rm,rb),(4,ins,im,ib)):
        for b in bb.values():
            rr=C.receiver(b['bbox'],sum(mass[i]for i in b['ids']),D)
            for k in (0,1):
                if rr[k]:layers[offset+k].append((rr[k],1))
    labelrows=[]
    for j,s in enumerate(labels):
        b=boxes[s];outside=[i for i in b['ids']if not S>>i&1];bm=sum((m[i]for i in outside),F(0));am=mm[s]-bm;core=C.receiver(b['bbox'],bm,D)
        if outside:
            tight=(min(p[i][0]for i in outside),max(p[i][0]for i in outside),min(p[i][1]for i in outside),max(p[i][1]for i in outside));tightcore=C.receiver(tight,bm,D)
        else:tightcore=[None,None]
        for k in (0,1):
            for offset,rr in ((6,full[s]),(8,core),(10,tightcore)):
                if rr[k]:layers[offset+k].append((rr[k],1<<j))
        labelrows.append(dict(mask=s,rank=s.bit_count(),mass=str(mm[s]),a=str(am),b=str(bm),full_bbox=list(map(str,b['bbox'])),full_span_squared=str(b['d2'])))
    b=boxes[S]['bbox'];q=[v*D for v in b];assert all(v.denominator==1 for v in q);a,z,yy,t=map(int,q)
    rr=A.G.root_floor(3*L.numerator*D2,4*L.denominator);rhos=(max(0,rr-1),rr+1);K=[(z-r,a+r,t-r,yy+r)for r in rhos]
    for k in (0,1):assert K[k][0]<K[k][1]and K[k][2]<K[k][3];layers[12+k].append((K[k],1))
    krecord=dict(formula='[xmax-rho,xmin+rho]x[ymax-rho,ymin+rho],rho=sqrt(3L)/2',inner=[str(F(v,D))for v in K[0]],outer=[str(F(v,D))for v in K[1]],rho_lower=str(F(rhos[0],D)),rho_upper=str(F(rhos[1],D)),rho_squared='3L/4');assert krecord==saved['K']
    values,zero,cells=area_grid(layers)
    def enc(lo,hi,signed=False):
        assert lo<=hi and(signed or lo>=0)
        return dict(lower=str(F(lo,D2)),upper=str(F(hi,D2)),normalized_lower=str(F(lo,D2)/L),normalized_upper=str(F(hi,D2)/L))
    parts=('global','central','peripheral');regions={}
    for name,v in values.items():
        regions[name]={part:enc(*v[j])for j,part in enumerate(parts)};g,k,o=v;res=(g[0]-k[1]-o[1],g[1]-k[0]-o[0]);assert res[0]<=0<=res[1];regions[name]['partition_identity_residual']=enc(*res,signed=True)
        assert regions[name]==saved['regions'][name],(data['id'],name)
    assert set(regions)==set(saved['regions'])
    identities={};net={}
    for j,part in enumerate(parts):
        dl,du=values['Delta'][j];il,iu=values['I'][j];cl,cu=values['internal_old_credit'][j];pl,pu=values['P_pure_mixed'][j]
        res=(dl-(iu-cl+pu),du-(il-cu+pl));assert res[0]<=0<=res[1];identities['Delta_minus_rhs_'+part]=enc(*res,signed=True)
        nlo,nhi=pl-cu,pu-cl;net[part]=dict(P_minus_credit=enc(nlo,nhi,signed=True),P_at_most_L_over_5_certified=F(pu,D2)<=L/5,P_strictly_above_L_over_5_certified=F(pl,D2)>L/5,net_at_most_L_over_5_certified=F(nhi,D2)<=L/5,net_strictly_above_L_over_5_certified=F(nlo,D2)>L/5)
        for state,surplus in (('W','W_old_surplus'),('Wprime','Wprime_old_surplus')):
            wl,wu=values[state][j];sl,su=values[surplus][j];res=(wl-du-su,wu-dl-sl);assert res[0]<=0<=res[1];identities[state+'_minus_rhs_'+part]=enc(*res,signed=True)
    assert identities==saved['area_identities']and net==saved['mixed_mass_gap_diagnostics']
    rcentral=values['R'][1];assert rcentral[0]==0 and (rcentral[1]==0)==saved['central_residual_halo_geometrically_certified_empty']
    lo=max(values[name][1][0]for name in ('Hfull','Delta','W','Wprime'));hi=min(values[name][1][1]for name in ('Hfull','Delta','W','Wprime'));assert lo<=hi and enc(lo,hi)==saved['central_Hfull_Delta_W_Wprime_common_interval']
    assert labelrows==sorted(saved['all_intersecting_labels'],key=lambda r:r['mask'])
    assert len(boxes)==saved['complete_closed_bbox_count']and len(elig)==saved['strict_eligible_count']and len(labels)==saved['eligible_intersect_S_count']and len(ib)==saved['internal_closed_bbox_count']and len(rb)==saved['residual_closed_bbox_count']
    assert all(F(v)==0 for v in saved['exact_zero_area_inclusion_checks'].values())
    comps=A.components(p,m,elig,boxes);assert sorted(c['mask']for c in comps)==sorted(saved['eligible_component_masks'])
    margin=A.G.span2(p,range(n))-F(35,33)*sum(m);assert margin==F(saved['whole_nonterminal_margin'])
    return dict(id=data['id'],S_mask=S,L=str(L),maximum_masks=winners,source_count=n,labels=len(labels),closed_bboxes=len(boxes),components=comps,K=krecord,regions=regions,mixed_mass_gap_diagnostics=net,zero_inclusion_areas=zero,compressed_grid_cells=cells,wall_seconds=time.perf_counter()-start)

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--archive-root',type=Path,default=Path(__file__).resolve().parents[2]);ap.add_argument('--spec',type=Path,default=Path(__file__).with_name('r069-independent-replay-inputs.json'));ap.add_argument('--output',type=Path,required=True);ap.add_argument('--workers',type=int,default=8);args=ap.parse_args()
    need(not sys.flags.optimize and Path('/Volumes/data').is_mount(),'Assertions and data volume required');need(1<=args.workers<=12,'workers1..12');root=args.archive_root.resolve()
    if root.name=='planar-centered-square-weak11':root=root.parent
    out=args.output.resolve();receipt=out.with_suffix('.receipt.json');need(out.is_relative_to(Path('/Volumes/data/CodexTemp').resolve())and out.suffix=='.json'and not out.exists()and not receipt.exists(),'New dataTemp output required');out.parent.mkdir(parents=True,exist_ok=True)
    start=time.perf_counter();meta=dict(pid=os.getpid(),python=sys.version,entry_sha256=sha(Path(__file__)),spec_sha256=sha(args.spec),archive_root=str(root),workers=args.workers,command=sys.argv)
    try:
        spec=json.loads(args.spec.read_text());pins=spec['dependencies']
        for rel,h in CORE_PINS.items():need(pins[rel]==h,'Frozen independent helper pin')
        for rel,h in pins.items():
            p=(root/rel).resolve();need(p.is_relative_to(root)and p.is_file()and sha(p)==h,f'Archive dependency mismatch {rel}')
        jobs=[]
        for c in spec['jobs']:
            need(c['input']in pins and c['output']in pins,'Unpinned job input/output');ip,op=root/c['input'],root/c['output'];data=json.loads(ip.read_text());saved=json.loads(op.read_text())
            need(saved['spec_sha256']=='5578573cedc914ec6ab57ea886bd2a16cd3646cdaa16a0205150f8383030535d','Frozen producer spec metadata')
            need(saved['status']=='PASS'and saved['source_sha256']==PRODUCER_SHA and saved['input_sha256']==sha(ip),'Production input/source chain');need(saved['geometry_source_sha256']=='4ef61b815b9c4a9a8972cdaf31716ae7917b59243a29fb8697c8042137788246'and saved['rectclip_source_sha256']=='4dc4fb494736c40334359f9950f052ee2e4a53b503883ffbc3f9bc0da69b1361','Frozen geometry metadata');jobs.append((data,saved))
        need(len(jobs)==8 and len({d['id']for d,r in jobs})==8,'Exactly8 distinct authorized states')
        if args.workers==1:boot(root);rows=[evaluate(j)for j in jobs]
        else:
            with ProcessPoolExecutor(max_workers=args.workers,initializer=boot,initargs=(str(root),))as pool:rows=list(pool.map(evaluate,jobs))
        result=dict(status='PASS_R69_INDEPENDENT_SOURCE_AND_CENTRAL_AREA_PARTITIONS',**meta,elapsed_seconds=time.perf_counter()-start,cases=rows,statistics=dict(states=len(rows),labels=sum(r['labels']for r in rows),compressed_grid_cells=sum(r['compressed_grid_cells']for r in rows)),dependency_pins=pins,scope='Exactly8 fixed states,11sets×global/central/peripheral; independent14layer2D difference/bitset with source reconstruction and directed K. No producer import; no new .2/general3.5 assertion.')
        out.write_text(json.dumps(result,indent=2)+'\n');meta.update(status=result['status'],elapsed_seconds=result['elapsed_seconds'],output_sha256=sha(out),statistics=result['statistics']);receipt.write_text(json.dumps(meta,indent=2)+'\n');print(json.dumps(meta,indent=2))
    except BaseException as e:
        meta.update(status='FAIL_R69_INDEPENDENT_REPLAY',error=str(e));receipt.write_text(json.dumps(meta,indent=2)+'\n');raise
if __name__=='__main__':main()
