#!/usr/bin/env python3
"""Independent rational-coordinate W/Wprime geometry; pinned archive reads only."""
from pathlib import Path
from fractions import Fraction as F
import argparse,hashlib,importlib.util,json,math,os,sys,time
sys.dont_write_bytecode=True
PREFIX='planar-centered-square-weak11/computation/'
DEPENDENCIES={PREFIX+'review_r067_erosion_states.py':'513e59dad17949c419ba98ec9c0e9ef3738871b92f32074df35eed874fd48eb0',PREFIX+'review_r066_grid_budgets.py':'c6524c436d8efd15a1b941c3e37448608622f6723fc51e3325a6eaa44aa159b1',PREFIX+'review_r064_mass_variation.py':'267b430e167801314bba81200250945bedcb1473918992f57ac03c1feb94cad4'}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def require(v,msg):
    if not v:raise RuntimeError(msg)
def boot(root):
    global B,A
    q=root/(PREFIX+'review_r067_erosion_states.py');s=importlib.util.spec_from_file_location('own_frozen_r067_grid',q);B=importlib.util.module_from_spec(s);s.loader.exec_module(B)
    B.boot(str(root/(PREFIX+'review_r066_grid_budgets.py')));A=B.A

def receiver(b,m,D):
    if m<=0 or max(b[1]-b[0],b[3]-b[2])**2>=m:return [None,None]
    k=A.G.root_floor(m.numerator*D*D,4*m.denominator)
    q=[z*D for z in b];assert all(z.denominator==1 for z in q);a,z,c,t=map(int,q)
    out=[]
    for r in (max(0,k-1),k+1):
        box=(z-r,a+r,t-r,c+r);out.append(box if box[0]<box[1]and box[2]<box[3]else None)
    return out

def five_source_certificate(p,m,S):
    eta=F(1,4096);expected=[(F(0),F(0)),(F(9,10),F(9,10)),(F(0),F(3,10)),(F(3,10),F(9,10)),(F(-1,2),F(0))]
    assert p==expected and m==[F(1,4)-2*eta,F(3,4),eta,eta,F(1,2)]and S==15
    # Rectangular patches, all endpoints rational and all comparisons exact.
    P=(F(-13,32),F(-25,64),F(-1,128),F(1,128))
    Q=(F(-7,32),F(-3,16),F(-1,64),F(1,64))
    def maxdist(point,box):return max(abs(point[0]-box[0]),abs(point[0]-box[1]),abs(point[1]-box[2]),abs(point[1]-box[3]))
    def mindist(point,box):
        dx=max(box[0]-point[0],F(0),point[0]-box[1]);dy=max(box[2]-point[1],F(0),point[1]-box[3]);return max(dx,dy)
    # N={A,C,q}: one common radius captures all three everywhere on P.
    N=(1<<0)|(1<<2)|(1<<4);r=F(27,64);mn=sum(m[i]for i in range(5)if N>>i&1)
    assert all(maxdist(p[i],P)<r for i in (0,2,4))and 4*r*r<mn
    # No eligible witness at P can ever capture B or D, even using all mass.
    assert all(4*mindist(p[i],P)**2>sum(m)for i in (1,3))
    # Full-N residual core excludes P because its A/C x=0 face is too far;
    # actual residual core is the q square, which contains P.
    assert 4*min(abs(P[0]),abs(P[1]))**2>m[4]
    assert 4*maxdist(p[4],P)**2<m[4]
    # No S-only label can meet P. Without B its mass <=1/4, while
    # all S positions have x>=0; with B the preceding total-mass test applies.
    assert 4*min(abs(P[0]),abs(P[1]))**2>sum(m[i]for i in (0,2,3))
    # Every eligible crossing label has B=q alone; thus all its Wprime
    # pieces avoid P. The internal A singleton creates old-halo surplus Q.
    assert 4*maxdist(p[0],Q)**2<m[0]and 4*maxdist(p[4],Q)**2<m[4]
    pa=(P[1]-P[0])*(P[3]-P[2]);qa=(Q[1]-Q[0])*(Q[3]-Q[2])
    assert pa==F(1,4096)and qa==F(1,1024)
    assert F(49,25)-F(35,33)*F(3,2)==F(203,550)>0
    return dict(strict_refinement_patch=list(map(str,P)),refinement_area_lower=str(pa),old_halo_surplus_patch=list(map(str,Q)),old_halo_surplus_area_lower=str(qa),actual_complete_source_mask=N,common_witness_radius=str(r),common_witness_area=str(4*r*r),complete_source_mass=str(mn),outside_sources_ruled_out_by_total_mass=True,nonterminal_margin='203/550')

def evaluate(data,bits):
    p=[tuple(map(F,v))for v in data['points']];m=list(map(F,data['initial_masses']));S=data['S_mask'];n=len(p)
    assert len(m)==n and min(m)>0 and len(set(p))==n and 0<S<(1<<n)-1
    D=math.lcm(1<<bits,*(v.denominator for pt in p for v in pt));D2=D*D
    boxes=A.G.boxes(p);mm={s:sum(m[i]for i in b['ids'])for s,b in boxes.items()}
    eligible=[s for s,b in boxes.items()if mm[s]>b['d2']];L=max(mm[s]for s in eligible);winners=sorted(s for s in eligible if mm[s]==L)
    assert S in winners and sum(m[i]for i in range(n)if S>>i&1)==L
    labels=sorted(s for s in eligible if s&S);full={s:receiver(boxes[s]['bbox'],mm[s],D)for s in eligible}
    keep=[i for i in range(n)if not S>>i&1];rp=[p[i]for i in keep];rm=[m[i]for i in keep];rb=A.G.boxes(rp)
    remain=[receiver(b['bbox'],sum(rm[i]for i in b['ids']),D)for b in rb.values()];layers=[[]for _ in range(12)];meta=[]
    for rr in full.values():
        for k in (0,1):
            if rr[k]:layers[k].append((rr[k],1))
    for rr in remain:
        for k in (0,1):
            if rr[k]:layers[2+k].append((rr[k],1))
    for j,s in enumerate(labels):
        b=boxes[s];outside=[i for i in b['ids']if not S>>i&1];bm=sum((m[i]for i in outside),F(0));am=mm[s]-bm
        core=receiver(b['bbox'],bm,D)
        if outside:
            tight=(min(p[i][0]for i in outside),max(p[i][0]for i in outside),min(p[i][1]for i in outside),max(p[i][1]for i in outside));tightcore=receiver(tight,bm,D)
        else:tightcore=[None,None]
        for k in (0,1):
            for offset,rr in ((4,full[s]),(6,core),(8,tightcore)):
                if rr[k]:layers[offset+k].append((rr[k],1<<j))
        meta.append(dict(mask=s,mass=str(mm[s]),a=str(am),b=str(bm)))
    for i in range(n):
        if S>>i&1:
            rr=receiver((p[i][0],p[i][0],p[i][1],p[i][1]),L,D)
            for k in (0,1):layers[10+k].append((rr[k],1))
    values,ind,zero,cells=B.grid(layers,len(labels))
    def interval(lo,hi):
        assert 0<=lo<=hi
        return dict(lower=str(F(lo,D2)),upper=str(F(hi,D2)),normalized_lower=str(F(lo,D2)/L),normalized_upper=str(F(hi,D2)/L))
    names=('W','Wprime','true_deletion','full_halo','remaining_halo','W_minus_Delta_overpayment','Wprime_minus_Delta_overpayment','W_minus_Wprime_refinement')
    metrics={name:interval(values[2*j],values[2*j+1])for j,name in enumerate(names)}
    for j,row in enumerate(meta):
        row['W_label_ring']=interval(ind[0][j],ind[1][j]);row['Wprime_label_ring']=interval(ind[2][j],ind[3][j]);assert F(ind[0][j],D2)<=F(row['a'])and F(ind[2][j],D2)<=F(row['a'])
    comps=A.components(p,m,eligible,boxes);b=boxes[S]['bbox'];fat=16*(b[1]-b[0])**2>9*L and 16*(b[3]-b[2])**2>9*L
    subsets=[]
    if n<=12:
        for s in range(1,1<<n):
            ids=[i for i in range(n)if s>>i&1];mass=sum(m[i]for i in ids);d2=A.G.span2(p,ids)
            subsets.append(dict(mask=s,mass=str(mass),span_squared=str(d2),strict_eligible=d2<mass))
        assert max(F(r['mass'])for r in subsets if r['strict_eligible'])==L
        assert sorted(r['mask']for r in subsets if r['strict_eligible']and F(r['mass'])==L)==winners
    certificate=None
    if data['id']=='r068-rank4-fat-strict-refinement-and-surplus':
        certificate=five_source_certificate(p,m,S)
        assert F(metrics['W_minus_Wprime_refinement']['lower'])>=F(1,4096)
        assert F(metrics['Wprime_minus_Delta_overpayment']['lower'])>=F(1,1024)
        assert winners==[15]and len(comps)==1 and fat and S.bit_count()==4
    return dict(id=data['id'],source_count=n,S_mask=S,S_rank=S.bit_count(),L=str(L),maximum_masks=winners,connected=len(comps)==1,whole_nonterminal=33*A.G.span2(p,range(n))>35*sum(m),both_axes_fat=fat,closed_bbox_count=len(boxes),eligible_closures=len(eligible),intersect_S_labels=len(labels),bits=bits,coordinate_area_denominator=str(D2),compressed_grid_cells=cells,metrics=metrics,zero_inclusion_areas=zero,labels=meta,all_nonempty_subsets=subsets,analytic_certificate=certificate)

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--archive-root',type=Path,default=Path(__file__).resolve().parents[2]);ap.add_argument('--spec',type=Path,default=Path(__file__).with_name('r068-independent-replay-inputs.json'));ap.add_argument('--output',type=Path,required=True);ap.add_argument('--bits',type=int,default=72);args=ap.parse_args()
    require(sys.flags.optimize==0 and args.bits>=16,'Assertions required; bits>=16');require(Path('/Volumes/data').is_mount(),'data volume required')
    root=args.archive_root.resolve()
    if root.name=='planar-centered-square-weak11':root=root.parent
    out=args.output.resolve();require(out.is_relative_to(Path('/Volumes/data/CodexTemp').resolve())and out.suffix=='.json','Output must be new dataTemp JSON');rp=out.with_suffix('.receipt.json');require(not out.exists()and not rp.exists(),'Use new output paths');out.parent.mkdir(parents=True,exist_ok=True)
    start=time.perf_counter();receipt=dict(pid=os.getpid(),python=sys.version,entry_sha256=sha(Path(__file__)),spec_sha256=sha(args.spec),archive_root=str(root),command=sys.argv)
    try:
        spec=json.loads(args.spec.read_text());pins=spec['dependencies']
        for rel,h in DEPENDENCIES.items():require(pins[rel]==h,'Frozen helper pin mismatch')
        for rel,h in pins.items():
            q=(root/rel).resolve();require(q.is_relative_to(root)and q.is_file()and sha(q)==h,f'Archive dependency mismatch {rel}')
        boot(root);rows=[]
        for c in spec['cases']:
            require(c['input']in pins,'Unpinned case');r=evaluate(json.loads((root/c['input']).read_text()),args.bits)
            if 'expected_result'in c:
                require(c['expected_result']in pins,'Unpinned expected result');expected=json.loads((root/c['expected_result']).read_text())
                for k,v in r['metrics'].items():require(v==expected[k],f'Expected interval mismatch {k}')
                r['producer_saved_result_comparison']='ALL_EIGHT_METRIC_INTERVALS_EXACT_MATCH'
            rows.append(r)
        result=dict(status='PASS_R68_INDEPENDENT_2D_STATE_GEOMETRY',**receipt,elapsed_seconds=time.perf_counter()-start,cases=rows,dependency_pins=pins,scope='Independent exact rational-coordinate12layer difference-grid/label-bitset states. Strict refinement and old-halo surplus do not refute either3.5 selection candidate or true deletion.')
        out.write_text(json.dumps(result,indent=2)+'\n');receipt.update(status=result['status'],elapsed_seconds=result['elapsed_seconds'],output_sha256=sha(out),cases=len(rows),labels=sum(r['intersect_S_labels']for r in rows));rp.write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))
    except BaseException as e:
        receipt.update(status='FAIL_R68_REPLAY',error=str(e));rp.write_text(json.dumps(receipt,indent=2)+'\n');raise
if __name__=='__main__':main()
