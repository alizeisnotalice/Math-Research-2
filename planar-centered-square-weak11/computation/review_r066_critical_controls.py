#!/usr/bin/env python3
"""Independent supplementary controls; imports only our independent grid audit."""
from fractions import Fraction as F
import hashlib,importlib.util,json,itertools,sys,time
from pathlib import Path
HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('own_r66_grid',HERE/'review_r066_grid_budgets.py')
A=importlib.util.module_from_spec(spec);spec.loader.exec_module(A);A.init_worker()

def all_subsets(case):
    pp=[tuple(map(F,p))for p in case['points']];mm=list(map(F,case['masses']));rows=[]
    for mask in range(1,1<<len(pp)):
        ids=[i for i in range(len(pp))if mask>>i&1]
        mass=sum(mm[i]for i in ids);span=max(max(pp[i][j]for i in ids)-min(pp[i][j]for i in ids)for j in (0,1))
        rows.append(dict(mask=mask,mass=str(mass),span_squared=str(span**2),strict_eligible=span**2<mass))
    return rows

def main():
    start=time.perf_counter();pins=[];rows=[]
    controls=[('THIN_CONTROL',A.EROSION/'five-source-control'),('FAT_CONTROL',A.EROSION/'three-source-fat-control')]
    for role,p in controls:
        bp,ep=p/'budget.json',p/'erosion.json';d=json.loads(bp.read_text());e=json.loads(ep.read_text())
        assert d['status']==e['status']=='PASS'and d['source_sha256']==A.BSHA and e['source_sha256']==A.ESHA
        assert e['budget_result_sha256']==A.sha(bp)
        c=d['cases'][0];r=A.evaluate_case((role,c,e['cases'][0]));r['all_nonempty_subsets']=all_subsets(c)
        r['independent_maximum_subset_masks']=[v['mask']for v in r['all_nonempty_subsets']if v['strict_eligible']and F(v['mass'])==1]
        assert r['independent_maximum_subset_masks']==r['all_maximum_masks'];rows.append(r)
        pins.extend([A.pin(bp),A.pin(ep)])
        for n in ('receipt.json','input.json','input-adapter.json'):
            if (p/n).exists():pins.append(A.pin(p/n))
    thin,fat=rows;tb,fb=thin['blocks'][0],fat['blocks'][0]
    assert not tb['fat_both_axes']and tb['strong_strict_fail']and tb['budget_pass']
    assert fat['connected']and fat['whole_nonterminal']and fat['all_maximum_masks']==[3]and fb['fat_both_axes']
    assert F(fb['budget']['lower'])>F(7,2)and F(fb['deletion']['upper'])<1 and fb['strong_strict_fail']
    # Entirely rational lower substitutions, separately checked by squaring.
    alpha,eps=F(1023,1024),F(1,1024)
    sl,tu=F(433,250),F(97,250)
    assert sl>0 and sl**2<3 and tu>0 and tu**2>F(3,20)
    blow=alpha**2+F(93,80)-eps+alpha*sl-tu-F(7,2)
    glow=F(33,80)-eps-tu-(F(7,4)-sl)
    assert blow==F(246397,131072000)>0 and glow==F(707,128000)>0
    assert (alpha+F(1,8))**2-F(35,33)*(1+eps)==F(6982433,34603008)>0
    # Five-source strict radical sign without evaluating a radical.
    assert 3*428**2-741**2==471>0
    # General erosion function checked at rational-square endpoints: u,v are
    # exact, sqrt(b),sqrt(m) chosen rational, hence every boundary value exact.
    scalar=0
    for u,v,rb,rm in itertools.product([F(0),F(1,4),F(1,2),F(1),F(3,2)],repeat=4):
        if rb>rm:continue
        def area(r):return max(F(0),r-u)*max(F(0),r-v)
        assert 0<=area(rm)-area(rb)<=rm**2-rb**2;scalar+=1
    # Independent generic bbox enumeration (no family-specific census shortcut)
    # at the smallest fat family input; no area or general-budget claim.
    N=32;pp=[(F(7*i,8*N),F(sign*7,16))for i in range(N)for sign in (-1,1)]+[(F(39,32),F(0))]
    mm=[F(1,2*N)]*(2*N)+[F(1,4)];boxes=A.G.boxes(pp)
    eligible=[(s,sum(mm[i]for i in b['ids']))for s,b in boxes.items()if sum(mm[i]for i in b['ids'])>b['d2']]
    best=max(m for s,m in eligible);winners=[s for s,m in eligible if m==best]
    assert best==1 and winners==[(1<<(2*N))-1]
    caps=[]
    for k in range(N//4+1,15*N//32):
        ids=list(range(2*k,2*N))+[2*N];mask=sum(1<<i for i in ids)
        assert mask in boxes and (mask,sum(mm[i]for i in ids))in eligible
        caps.append(F(N-k,N))
    assert sum(caps)==F(41,64)*(7*N//32-1)
    out=dict(status='PASS_INDEPENDENT_R66_CRITICAL_CONTROLS',script_sha256=A.sha(Path(__file__)),python=sys.version,
             elapsed_seconds=time.perf_counter()-start,source_pins=pins+[A.pin(HERE/'review_r066_grid_budgets.py'),A.pin(HERE/'review_r064_mass_variation.py')],
             cases=rows,rational_radical_certificate=dict(B_excess_lower=str(blow),strong_excess_lower=str(glow),sqrt3_lower=str(sl),sqrt_3over20_upper=str(tu)),
             scalar_erosion_cases=scalar,fat_repetition_N32=dict(bbox_count=len(boxes),eligible_count=len(eligible),maximum_mass=str(best),maximizers=len(winners),selected_tails=len(caps),capacity_sum=str(sum(caps))),
             scope='Separate controls; original 201-block replay unchanged. New fat example refutes fixed theta 3/4 envelope selection, not true deletion selection or C2.')
    p=HERE/'results/round066-independent-critical-controls.json';p.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items()if k not in ('source_pins','cases')},indent=2))
if __name__=='__main__':main()
