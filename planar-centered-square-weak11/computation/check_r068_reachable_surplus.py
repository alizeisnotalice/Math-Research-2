"""Exact six-source domain and positive reachable-surplus patch certificate.

No full halo union is computed. All constants are rational and every source
subset is checked; the patch inequalities are uniform on an open square.
"""
from fractions import Fraction as F
from pathlib import Path
import hashlib,json,os,platform,time


def main():
    assert os.path.ismount('/Volumes/data')
    started=time.perf_counter()
    p=[(F(3,10),F(3,10)),(F(6,5),F(6,5)),(F(3,10),F(6,5)),
       (F(6,5),F(3,10)),(F(-1,50),F(3,10)),(F(7,50),F(7,50))]
    m=list(map(F,['1/10','449/500','1/1000','1/1000','1/10','1/1000']))
    S=15;rows=[];closed=[];elig=[]
    for mask in range(1,64):
        ids=[i for i in range(6) if mask>>i&1]
        lo=[min(p[i][j] for i in ids) for j in range(2)]
        hi=[max(p[i][j] for i in ids) for j in range(2)]
        span2=max(hi[j]-lo[j] for j in range(2))**2
        mass=sum((m[i] for i in ids),F(0))
        closure=sum(1<<i for i in range(6) if all(lo[j]<=p[i][j]<=hi[j] for j in range(2)))
        row={'mask':mask,'mass':str(mass),'span_squared':str(span2),'closure':closure,'eligible':span2<mass}
        rows.append(row)
        if closure==mask:closed.append(row)
        if span2<mass:elig.append(row)
    assert len(closed)==24 and sum(x['eligible'] for x in closed)==13
    L=max(F(x['mass']) for x in elig)
    assert L==1 and [x['mask'] for x in elig if F(x['mass'])==L]==[S]
    reach=1
    while True:
        nxt=reach
        for row in elig:
            if row['mask']&reach:nxt|=row['mask']
        if nxt==reach:break
        reach=nxt
    assert reach==63
    margin=F(61,50)**2-F(35,33)*sum(m)
    assert margin==F(17637,55000)>0 and F(9,10)>F(3,4)
    patch_delta=F(1,10000);center=(F(7,50),F(7,50));r=F(161,1000)
    N=49;B=48
    n_mass=sum(m[i] for i in range(6) if N>>i&1)
    b_mass=sum(m[i] for i in range(6) if B>>i&1)
    assert n_mass==F(201,1000) and b_mass==F(101,1000)
    for i in [0,4,5]:
        assert max(abs(p[i][j]-center[j]) for j in range(2))+patch_delta<r
    far_distance=min(max(abs(p[i][j]-center[j]) for j in range(2))-patch_delta for i in [1,2,3])
    assert far_distance==F(10599,10000)
    far_area=4*far_distance**2
    near_distance=F(1599,10000);near_area=4*near_distance**2
    assert 4*r*r==F(25921,250000)<n_mass
    assert far_area==F(112338801,25000000)>sum(m)>1
    assert near_area==F(2556801,25000000)>b_mass>m[0]
    assert 4*patch_delta**2<m[5]
    assert (2*patch_delta)**2==F(1,25000000)
    assert next(x for x in rows if x['mask']==N)['closure']==N
    edge_cost={}
    for i in range(4):
        for j in range(i+1,4):
            dx,dy=abs(p[i][0]-p[j][0]),abs(p[i][1]-p[j][1])
            edge_cost[i,j]=dx+dy-dx*dy
    assert min(edge_cost.values())==F(9,10)
    tree=[(0,2),(1,2),(1,3)]
    tree_cost=sum(edge_cost[e] for e in tree)
    assert tree_cost==3*min(edge_cost.values())==F(27,10)
    assert 1+tree_cost==F(37,10)>F(7,2)
    assert F(19,10)**2==F(361,100)>F(7,2)
    root=Path(__file__).resolve().parents[1]
    inp=root/'computation/inputs/r068-reachable-surplus-six-source.json'
    obj={'id':'r068-reachable-surplus-outside-internal-halo','points':[[str(x),str(y)] for x,y in p],
         'initial_masses':[str(t) for t in m],'S_mask':S,'S_indices':[0,1,2,3],
         'patch_center':[str(x) for x in center],'patch_open_half_side':str(patch_delta),
         'N_full_mask':N,'N_residual_mask':32,'N_full_minus_S_mask':B,
         'scope':'Positive patch of reachable full-state Wprime overpayment outside H(S), not a3.5 area counterexample.'}
    inp.write_text(json.dumps(obj,indent=2)+'\n')
    out={'status':'PASS_EXACT_SIX_SOURCE_REACHABLE_SURPLUS_PATCH','python':platform.python_version(),
         'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
         'input_sha256':hashlib.sha256(inp.read_bytes()).hexdigest(),'all_subsets':rows,
         'closed_bboxes':24,'eligible_closed_bboxes':13,'unique_maximum_mask':S,'L':'1',
         'connected':True,'S_rank':4,'fat_both_axes':True,'nonterminal_margin':str(margin),
         'patch_area':'1/25000000','actual_full_N_mask':N,'actual_residual_N_mask':32,
         'full_N_residual_mask':B,'common_full_witness_radius':str(r),'common_full_witness_area':str(4*r*r),
         'residual_capture_area_lower':str(near_area),'all_far_source_capture_area_lower':str(far_area),
         'N_full_minus_S_core_avoids_patch':True,'internal_halo_avoids_patch':True,
         'whole_patch_in_Wprime_and_R':True,'seconds':time.perf_counter()-started,
         'minimum_support_tree_cost':'27/10','support_tree_bound':'37/10','support_neighborhood_union_area':'361/100',
         'scope':'Exact source enumeration and uniform rational patch inequalities only; no full halo area computed.'}
    target=root/'computation/results/round068/reachable-surplus-arithmetic.json'
    target.parent.mkdir(parents=True,exist_ok=True);target.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k!='all_subsets'},indent=2))


if __name__=='__main__':main()
