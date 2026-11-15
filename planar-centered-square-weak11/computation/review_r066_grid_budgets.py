#!/usr/bin/env python3
"""Independent R66 10-layer 2D difference-grid budget/erosion replay.

No producer imports or strip/interval-merge/subtract-rectangle algorithms.
Per-high-label bitsets compute unions of rectangle DIFFERENCES directly.
"""
from collections import Counter,defaultdict
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction as F
import argparse,hashlib,importlib.util,json,os
from pathlib import Path
import sys,time

HERE=Path(__file__).resolve().parent
ROOT=Path('/Volumes/data/CodexProjects/T2')
FRONT=ROOT/'.local/compute/frontier066'
COMPUTE=ROOT/'.local/workspaces/pcw-position-slab-impl/planar-centered-square-weak11'
EROSION=COMPUTE/'computation/results/round066-erosion'
BSHA='211d33e21e47e6e14290e4ddbc3b3a36e247131de1bc53114091b7766d0fe8b6'
ESHA='4dc4fb494736c40334359f9950f052ee2e4a53b503883ffbc3f9bc0da69b1361'

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def pin(p):return dict(path=str(p),sha256=sha(p),bytes=p.stat().st_size)


def init_worker():
    global G
    spec=importlib.util.spec_from_file_location('own_r64_independent_geometry',HERE/'review_r064_mass_variation.py')
    G=importlib.util.module_from_spec(spec);spec.loader.exec_module(G)


def rect(bounds,mass,bits):
    if mass<=0 or max(bounds[1]-bounds[0],bounds[3]-bounds[2])**2>=mass:return [None,None]
    D=1<<bits;k=G.root_floor(mass.numerator*D*D,4*mass.denominator)
    q=[v*D for v in bounds];assert all(v.denominator==1 for v in q)
    a,b,c,d=map(int,q);out=[]
    for r in (max(0,k-1),k+1):
        z=(b-r,a+r,d-r,c+r);out.append(z if z[0]<z[1]and z[2]<z[3]else None)
    return out


def components(points,masses,eligible,boxes):
    parent=list(range(len(points)))
    def find(i):
        while parent[i]!=i:parent[i]=parent[parent[i]];i=parent[i]
        return i
    for s in eligible:
        ids=boxes[s]['ids'];a=find(ids[0])
        for i in ids[1:]:parent[find(i)]=a
    groups=defaultdict(list)
    for i in range(len(points)):groups[find(i)].append(i)
    result=[]
    for ids in groups.values():
        mask=sum(1<<i for i in ids);m=sum(masses[i]for i in ids);d2=G.span2(points,ids)
        result.append(dict(mask=mask,rank=len(ids),mass=str(m),span_squared_over_mass=str(d2/m),terminal_35_33=33*d2<=35*m))
    return sorted(result,key=lambda r:r['mask'])


def grid(layers,nhigh):
    # Each item is(rectangle,integer weight). High-related weights are distinct
    # powers of2, one rectangle per label, hence prefix sums are exact bitsets.
    xs=sorted({r[j]for layer in layers for r,w in layer for j in (0,1)})
    ys=sorted({r[j]for layer in layers for r,w in layer for j in (2,3)})
    nx,ny=len(xs),len(ys);xi={v:i for i,v in enumerate(xs)};yi={v:i for i,v in enumerate(ys)}
    assert nx*ny<4_000_000,'independent grid memory guard'
    aa=[[0]*(nx*ny)for _ in layers]
    for ar,layer in zip(aa,layers):
        for r,w in layer:
            xl,xr,yl,yr=xi[r[0]],xi[r[1]],yi[r[2]],yi[r[3]]
            for x,y,v in ((xl,yl,w),(xr,yl,-w),(xl,yr,-w),(xr,yr,w)):ar[x*ny+y]+=v
    # metrics: oldi/o, highi/o, rawhi/o, holei/o, deltai/o, ringoutsidei/o
    totals=[0]*12;individual=[[0]*nhigh for _ in range(8)]
    for ix in range(nx):
        for iy in range(ny):
            j=ix*ny+iy
            for ar in aa:
                if ix:ar[j]+=ar[j-ny]
                if iy:ar[j]+=ar[j-1]
                if ix and iy:ar[j]-=ar[j-ny-1]
                assert ar[j]>=0
            if ix+1==nx or iy+1==ny:continue
            fi,fo,ri,ro,ci,co,mi,mo,bi,bo=[a[j]for a in aa]
            assert not fi or fo;assert not ri or ro;assert not ci or co
            assert not mi or fi;assert not mo or fo
            assert not bi&~mi and not bo&~mo
            assert not bi or ri;assert not bo or ro
            ringi,ring_o=mi&~bo,mo&~bi
            area=(xs[ix+1]-xs[ix])*(ys[iy+1]-ys[iy])
            flags=(ci and ri,co and ro,mi and not co and not ro,mo and not ci and not ri,
                   mi and not co,mo and not ci,ci and not fo,co and not fi,fi and not ro,fo and not ri,
                   ringi and not co,ring_o and not ci)
            for k,v in enumerate(flags):
                if v:totals[k]+=area
            masks=(ringi,ring_o,ringi if not co else 0,ring_o if not ci else 0,
                   mi if not co else 0,mo if not ci else 0,
                   mi if not co and not ro else 0,mo if not ci and not ri else 0)
            for ar,bits in zip(individual,masks):
                while bits:
                    bit=bits&-bits;ar[bit.bit_length()-1]+=area;bits^=bit
    return totals,individual,(nx-1)*(ny-1),nx,ny


def evaluate_case(job):
    role,case,erosion=job;started=time.perf_counter()
    points=[tuple(map(F,p))for p in case['points']];masses=list(map(F,case['masses']));n=len(points)
    assert len(set(points))==n and min(masses)>0
    boxes=G.boxes(points);mm={s:sum(masses[i]for i in b['ids'])for s,b in boxes.items()}
    eligible=[s for s,b in boxes.items()if mm[s]>b['d2']];Lmax=max(mm[s]for s in eligible)
    winners=sorted(s for s in eligible if mm[s]==Lmax)
    assert winners==sorted(case['all_maximum_masks'])
    assert sorted(b['S_mask']for b in case['blocks'])==[s for s in winners if s!=(1<<n)-1]
    comps=components(points,masses,eligible,boxes);connected=len(comps)==1
    wholewide=33*G.span2(points,range(n))>35*sum(masses)
    emap={b['S_mask']:b for b in erosion['blocks']};out=[]
    for saved in case['blocks']:
        S=saved['S_mask'];L=mm[S];bits=saved['bits'];D=1<<bits;D2=D*D;theta=F(saved['theta'])
        assert L==Lmax==F(saved['L'])and 0<theta<1
        assert saved['components']==comps and saved['closed_bbox_count']==len(boxes)and saved['strict_eligible_count']==len(eligible)
        hi=sorted(s for s in eligible if s&S and mm[s]>theta*L);assert S in hi
        assert len(hi)==saved['high_competitor_count_including_S']
        full={s:rect(boxes[s]['bbox'],mm[s],bits)for s in eligible}
        keep=[i for i in range(n)if not S&(1<<i)];pp=[points[i]for i in keep];rm=[masses[i]for i in keep];rb=G.boxes(pp)
        remaining=[rect(b['bbox'],sum(rm[i]for i in b['ids']),bits)for b in rb.values()]
        a,b,c,d=boxes[S]['bbox'];rad=G.root_floor((theta*L).numerator*D*D,4*(theta*L).denominator)
        coords=[q*D for q in (a,b,c,d)];assert all(q.denominator==1 for q in coords);aD,bD,cD,dD=map(int,coords)
        C=[(aD-r,bD+r,cD-r,dD+r)for r in (max(0,rad-1),rad+1)]
        carea=[(r[1]-r[0])*(r[3]-r[2])for r in C]
        layers=[[]for _ in range(10)]
        for rr in full.values():
            for k in (0,1):
                if rr[k]:layers[k].append((rr[k],1))
        for rr in remaining:
            for k in (0,1):
                if rr[k]:layers[2+k].append((rr[k],1))
        layers[4]=[(C[0],1)];layers[5]=[(C[1],1)]
        ab=[];bm=[]
        for j,s in enumerate(hi):
            aval=sum(masses[i]for i in boxes[s]['ids']if S&(1<<i));bval=mm[s]-aval;ab.append((aval,bval))
            rr=rect(boxes[s]['bbox'],bval,bits);bm.append(rr)
            for k in (0,1):
                if full[s][k]:layers[6+k].append((full[s][k],1<<j))
                if rr[k]:layers[8+k].append((rr[k],1<<j))
        vals,ind,cells,nx,ny=grid(layers,len(hi))
        oldi,oldo,highi,higho,rawi,rawo,holei,holeo,deli,delo,ringi,ringo=vals
        bud_i=max(0,carea[0]-oldo+highi);bud_o=carea[1]-oldi+higho
        def interval(lo,up):
            assert 0<=lo<=up
            return dict(lower=str(F(lo,D2)),upper=str(F(up,D2)),normalized_lower=str(F(lo,D2)/L),normalized_upper=str(F(up,D2)/L))
        computed={name:interval(lo,up)for name,lo,up in (
            ('container_area',*carea),('low_old_credit',oldi,oldo),('high_unpaid_union',highi,higho),
            ('high_union_outside_container_before_old_credit',rawi,rawo),('low_hole',holei,holeo),
            ('true_deletion',deli,delo),('analytic_envelope_budget',bud_i,bud_o))}
        for k,v in computed.items():assert saved[k]==v,(case['id'],S,k)
        resid=(F(bud_i-delo-holeo,D2),F(bud_o-deli-holei,D2));assert resid[0]<=0<=resid[1]
        assert saved['identity_budget_minus_delta_minus_hole']==dict(lower=str(resid[0]),upper=str(resid[1]),contains_zero=True)
        score=(F(highi-oldo,D2)/L,F(higho-oldi,D2)/L)
        assert saved['high_unpaid_minus_low_old_credit_normalized']==dict(lower=str(score[0]),upper=str(score[1]))
        sq3=G.root_floor(3*D2,1);threshold=(F(7,4)-F(sq3+1,D),F(7,4)-F(sq3,D))
        assert saved['universal_compensation_threshold']==dict(lower=str(threshold[0]),upper=str(threshold[1]))
        strong=score[1]<=threshold[0];fail=score[0]>threshold[1];budgetpass=F(bud_o,D2)<=F(7,2)*L
        assert saved['budget_certifies_7over2']==budgetpass and saved['certifies_universal_compensation_condition']==strong
        ers=emap[S];assert ers['ring_union_outside_container']==interval(ringi,ringo)
        assert ers['actual_high_unpaid_union_from_budget']==computed['high_unpaid_union']
        ringrows={r['mask']:r for r in ers['competitors']};individual={r['mask']:r for r in saved['individual_competitors_DIAGNOSTIC_ONLY_NOT_UNION_SUM']}
        ring_over=0;individuals=[]
        for j,s in enumerate(hi):
            aval,bval=ab[j];rr=ringrows[s];ir=individual[s]
            assert (F(rr['a_intersection_mass']),F(rr['b_remaining_mass']),F(rr['total_mass']))==(aval,bval,mm[s])
            assert rr['erosion_ring_area']==interval(ind[0][j],ind[1][j])
            assert rr['erosion_ring_outside_container']==interval(ind[2][j],ind[3][j])
            assert F(ind[0][j],D2)<=aval
            assert rr['geometric_ring_upper_at_most_a']==(F(ind[1][j],D2)<=aval)
            area=[]
            for r in full[s]:area.append(0 if not r else (r[1]-r[0])*(r[3]-r[2]))
            assert ir['receiver_area']==interval(*area)
            assert ir['outside_container_area']==interval(ind[4][j],ind[5][j])
            assert ir['individually_unpaid_after_remaining_halo']==interval(ind[6][j],ind[7][j])
            ring_over+=F(ind[1][j],D2)>aval
            individuals.append(dict(mask=s,a=str(aval),ring_lower=str(F(ind[0][j],D2)),ring_upper=str(F(ind[1][j],D2))))
        fat=16*(b-a)**2>9*L and 16*(d-c)**2>9*L
        out.append(dict(S_mask=S,L=str(L),fat_both_axes=fat,strong_pass=strong,strong_strict_fail=fail,
                        budget_pass=budgetpass,score=list(map(str,score)),budget=computed['analytic_envelope_budget'],
                        deletion=computed['true_deletion'],hole=computed['low_hole'],ring_union=interval(ringi,ringo),
                        high_labels=len(hi),all_individual_ring_rows=individuals,ring_upper_rounding_over_a=ring_over,
                        compressed_grid_cells=cells,grid_shape=[nx,ny],S_is_whole_component=any(v['mask']==S for v in comps)))
    return dict(role=role,id=case['id'],source_count=n,all_maximum_masks=winners,proper_blocks=len(out),
                connected=connected,whole_nonterminal=wholewide,components=comps,blocks=out,wall_seconds=time.perf_counter()-started)


def jobs_and_pins():
    specs=[('R64',FRONT/'r64-final-theta.json',EROSION/'r64-final-erosion.json')]
    specs += [('R63',FRONT/f'stratum-{i:02}.json',EROSION/f'stratum-{i:02}-erosion.json')for i in range(12)]
    specs += [('THIN_CONTROL',EROSION/'five-source-control/budget.json',EROSION/'five-source-control/erosion.json')]
    jobs=[];pins=[]
    for role,bp,ep in specs:
        d=json.loads(bp.read_text());e=json.loads(ep.read_text())
        assert d['status']==e['status']=='PASS'and d['source_sha256']==BSHA and e['source_sha256']==ESHA
        assert e['budget_result_sha256']==sha(bp)
        emap={c['id']:c for c in e['cases']}
        for c in d['cases']:jobs.append((role,c,emap[c['id']]))
        pins += [pin(bp),pin(ep)]
    return jobs,pins


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True);ap.add_argument('--workers',type=int,default=12);ap.add_argument('--case-id')
    args=ap.parse_args();assert os.path.ismount('/Volumes/data')and 1<=args.workers<=12
    start=time.perf_counter();jobs,pins=jobs_and_pins()
    if args.case_id:jobs=[j for j in jobs if j[1]['id']==args.case_id];assert jobs
    if args.workers==1:
        init_worker();rows=[evaluate_case(j)for j in jobs]
    else:
        with ProcessPoolExecutor(max_workers=args.workers,initializer=init_worker)as pool:rows=list(pool.map(evaluate_case,jobs))
    blocks=[b for r in rows for b in r['blocks']];r63=[b for r in rows if r['role']=='R63'for b in r['blocks']]
    stats=dict(cases=len(rows),blocks=len(blocks),R63_cases=sum(r['role']=='R63'for r in rows),R63_blocks=len(r63),
               R63_budget_pass=sum(b['budget_pass']for b in r63),R63_strong_pass=sum(b['strong_pass']for b in r63),
               R63_strong_strict_fail=sum(b['strong_strict_fail']for b in r63),individual_high_labels=sum(b['high_labels']for b in blocks),
               individual_ring_upper_rounding_exceeds_a=sum(b['ring_upper_rounding_over_a']for b in blocks),
               total_compressed_grid_cells=sum(b['compressed_grid_cells']for b in blocks),max_compressed_grid_cells=max((b['compressed_grid_cells']for b in blocks),default=0))
    if not args.case_id:assert (stats['cases'],stats['blocks'],stats['R63_blocks'],stats['R63_strong_pass'],stats['R63_strong_strict_fail'])==(74,201,199,181,18)
    out=dict(status='PASS_INDEPENDENT_2D_GRID_R66_BUDGET_AND_RING_REPLAY',python=sys.version,executable=sys.executable,
             script_sha256=sha(Path(__file__)),elapsed_seconds=time.perf_counter()-start,workers=args.workers,statistics=stats,
             arithmetic='Exact integer/Fraction; own coordinate-tight bboxes, binary-search root bounds;10layer compressed2D difference-grid andper-label bitsets, no producer import/strip/merge/subtractrectangle.',
             source_pins=pins+[pin(HERE/'review_r064_mass_variation.py')],cases=rows,
             scope='Exact specified-input budget/ring intervals, all recorded proper maximum ties andclassification. Finite evidence only, not universal compensation/selection/deletion/C2 theorem.')
    args.out.parent.mkdir(parents=True,exist_ok=True);args.out.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items()if k not in ('source_pins','cases')},indent=2))


if __name__=='__main__':main()
