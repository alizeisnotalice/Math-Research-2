"""python3 -B check_interface_two.py. Uses adjacent check_joint_cut.py."""
from collections import Counter
from check_joint_cut import aff, cut_sparse, parity, rank, spectrum, walsh
import json
import platform


def audit():
    n=4
    size=2**n
    small={t:spectrum(t,2) for t in range(16)}
    stats=Counter()
    for t in range(2**size):
        active=[set(small[(t>>(4*i))&15])-{0} for i in range(4)]
        union=sorted(set().union(*active))
        if len(union)!=2 or any(len(a)>1 for a in active):
            continue
        signed={v:a for v,a in enumerate(walsh(1-2*((t>>x)&1) for x in range(size))) if a}
        sides=[{v for v in signed if (v&3)==b} for b in union]
        pure=[v for v in signed if v and (v&3)==0]
        d=rank(signed)
        if not any(rank(sides[0]|{p for j,p in enumerate(pure) if mask>>j&1})+
                   rank(sides[1]|{p for j,p in enumerate(pure) if not(mask>>j&1)})-d==2
                   for mask in range(2**len(pure))):
            continue
        a,b=({i for i,g in enumerate(active) if v in g} for v in union)
        key='-'.join(map(str,sorted((len(a),len(b)))))
        stats[key]+=1
        valid=[q for q in range(1,4) if
               {parity(q&i) for i in a}.isdisjoint({parity(q&i) for i in b})]
        if max(len(a),len(b))>2:
            assert not valid and sorted((len(a),len(b)))==[1,3]
            continue
        assert valid
        for q in valid:
            slices=cut_sparse(signed,(q<<2,))
            assert len(set(slices[0]) & set(slices[1]))<=2
            assert all(aff(S)<aff(signed) for S in slices)
            assert any((aff(S)+1)**2*len(signed)>=(aff(signed)+1)**2*len(S) for S in slices)
            stats['separating_single_cuts_checked']+=1
    return dict(stats)


def equal_pair_fixture():
    values=[]
    for x in range(64):
        address=x>>4
        values.append(1-2*((x>>address)&1))
    signed={v:a for v,a in enumerate(walsh(values)) if a}
    A={v for v in signed if v&15 in (1,2)}
    B=set(signed)-A
    assert rank(A)+rank(B)-rank(signed)==2
    rows=cut_sparse(signed,(32,))
    assert len(set(rows[0]) & set(rows[1]))==0
    assert all((aff(row)+1)**2*len(signed)>=(aff(signed)+1)**2*len(row) for row in rows)
    return {"variables":6,"sparsity":len(signed),"child_sparsities":list(map(len,rows)),
            "interface_dimension":2,"active_counts":[2,2]}


if __name__=='__main__':
    print(json.dumps({"python":platform.python_version(),"tolerance":0,"seed":None,
                      "scope":"all 65536 four-variable tables; fixed two address coordinates",
                      "counts":audit(),"equal_pair_example":equal_pair_fixture()},indent=2))
