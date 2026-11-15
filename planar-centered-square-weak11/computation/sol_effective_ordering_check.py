from fractions import Fraction as F
from pathlib import Path
import json,platform,hashlib,math
from itertools import combinations,permutations
nested=[]
for N in (3,8,32):
    delta=F(1,100*N*N);D=(N-1)*delta
    assert D*D<1 and D<F(1,4*N)
    nested.append(dict(N=N,delta=str(delta),diameter=str(D),lower_radius_gap=str(F(1,4*N)),full_label_load=N,envelope_area=N))
pointwise=[]
for delta in (F(1,10),F(1,100),F(1,1000)):
    eta=delta**4;d=1-delta;near=F(1,2)-delta;far=F(3,2)-2*delta
    assert d*d<1 and near>0
    small=eta/(4*near*near);allmass=1/(4*far*far)
    assert small<allmass
    pointwise.append(dict(delta=str(delta),eta=str(eta),bbox_span_squared=str(d*d),light_average=str(small),full_average=str(allmass),maximal_average=str(allmass)))
d=F(7,10);a=F(49,100)
assert (4*d-1)**2>2 and 4*d-1>0 and d*d<F(1,2)
assert (8*a-2)**2>2 and 8*a-2>0
points=[(x,y) for x in (-d/2,d/2) for y in (-a,a)]
def span(ids):
    pts=[points[i] for i in ids]
    return max(max(z[q] for z in pts)-min(z[q] for z in pts) for q in (0,1))
classifications=[]
for k in range(1,5):
    for S in combinations(range(4),k):
        valid=span(S)**2<F(k,4)
        expected=k in (1,4) or (k==2 and points[S[0]][1]==points[S[1]][1])
        assert valid==expected
        classifications.append(dict(indices=list(S),mass=str(F(k,4)),span_squared=str(span(S)**2),eligible=valid))
orders=[]
for order in permutations(range(4)):
    old=[];masses=[]
    for i in order:
        eff=F(0)
        for k in range(len(old)+1):
            for S in combinations(old,k):
                if span(S+(i,))**2<F(k+1,4):eff=max(eff,F(k,4))
        masses.append(eff+F(1,4));old.append(i)
    assert masses[-1]==1
    otherrow=[masses[j] for j in range(3) if points[order[j]][1]!=points[order[-1]][1]]
    lastrow=[masses[j] for j in range(3) if points[order[j]][1]==points[order[-1]][1]]
    assert otherrow==[F(1,4),F(1,2)] and lastrow==[F(1,4)]
    value=F(13,8)+d/2;assert value==F(79,40)
    orders.append(dict(order=list(order),proxy_masses=[str(z) for z in masses],union_area_by_proved_disjointness=str(value)))
proof=Path(__file__).resolve().parents[1]/'proofs/sol-effective-ordering-draft.md' 
print(json.dumps(dict(python=platform.python_version(),method='deterministic exact Fraction; single process, no random seed or tolerance',nested_examples=nested,pointwise_examples=pointwise,inherited_old_bound_times_9=str(9*F('3.590865312933')),rectangle_four_point=dict(d=str(d),subsets=classifications,orders=orders,limiting_proxy_formula='13/8+sqrt(2)/4',limiting_proxy_float_display=13/8+math.sqrt(2)/4,true_halo_formula='2-sqrt(2)*d+(1-d)/50',true_halo_float_display=2-math.sqrt(2)*float(d)+(1-float(d))/50),proof_sha256=hashlib.sha256(proof.read_bytes()).hexdigest()),indent=2))
