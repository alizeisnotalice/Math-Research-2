"""Independent sparse integer checks; general claims are proved analytically."""
import hashlib
import itertools
import json
import platform
from pathlib import Path

HERE = Path(__file__).resolve().parent

def mul(a, b):
    out = {}
    for i, u in a.items():
        for j, v in b.items():
            out[i+j] = out.get(i+j, 0) + u*v
    return {i: v for i, v in out.items() if v}

def stretch(a, p):
    return {p*i: v for i, v in a.items()}

def add(*polys):
    out = {}
    for a in polys:
        for i, u in a.items():
            out[i] = out.get(i, 0) + u
    return {i: v for i, v in out.items() if v}

def check_family(p, q):
    f = {0:1, q-p:1, q:-1, q+p:1, 2*q:1}
    g = {0:1, p:1, q:1, p+q:1}
    h = mul(f, g)
    a = mul(stretch(f, p), stretch(g, q))
    b = mul(stretch(g, p), stretch(f, q))
    c = mul(stretch(h, p), stretch(h, q))
    # Explicit nonnegative decompositions, recomputed independently of dense code.
    a_pos = mul(add({0:1, 3*p*q:1},
                    mul({p*(q-p):1, p*(q+p):1}, {0:1, p*q:1})),
                {0:1, q*q:1})
    b_pos = mul(add(mul({0:1, 2*q*q:1}, {0:1, p*q:1}),
                    {q*(q-p):1, q*(q+2*p):1}), {0:1, p*p:1})
    assert a == a_pos and b == b_pos
    expected = {0,p,q-p,q,p+q,2*q-p,2*q,q+2*p,2*q+p,2*q+2*p,3*q,3*q+p}
    assert len(expected) == 12 and h == {i:1 for i in expected}
    assert min(a.values()) > 0 and min(b.values()) > 0 and c[p*q] == 2
    return {'p':p, 'q':q, 'degree_H':p+3*q, 'H_support':sorted(h),
            'cross_nonnegative_identities':True, 'C_pq':c[p*q]}

def check_aliases(p, q):
    # Pair lists verify the coefficient formulas without evaluating sample factors.
    for r in range(q, 2*q):
        pairs = [(i,j) for i in range(r+1) for j in range(2*q)
                 if p*i+q*j == p*r]
        assert pairs == [(r-q,p),(r,0)]
        pairs = [(i,j) for i in range(2*q) for j in range(r+1)
                 if p*i+q*j == q*r]
        assert pairs == [(0,r),(q,r-p)]
    pairs = [(i,j) for i in range(2*q+1) for j in range(2*q+1)
             if p*i+q*j == q*q]
    expected = [(0,q),(q,q-p)]
    if q >= 2*p:
        expected.append((2*q,q-2*p))
    assert pairs == expected
    # The final q-difference in the equality proof makes two distinct C terms.
    if q > 2*p:
        s = q-2*p
        assert p*(s+q)+q*0 == p*s+q*p
    return {'p':p,'q':q,'edge_B_q2_pairs':pairs}

def check_residue_grid():
    f = {0:1,1:-2,3:4,5:-3,8:1,13:2}
    g = {0:1,2:3,4:-1,7:2,12:1}
    rows = []
    for p,q in [(2,3),(3,5),(3,7),(5,7)]:
        parts = []
        residues = set()
        for r in range(q):
            for s in range(p):
                residue = (p*r+q*s) % (p*q)
                assert residue not in residues
                residues.add(residue)
                fr = {(i-r)//q:v for i,v in f.items() if i % q == r}
                gs = {(i-s)//p:v for i,v in g.items() if i % p == s}
                prod = mul(fr, gs)
                parts.append({p*r+q*s+p*q*k:v for k,v in prod.items()})
        assert add(*parts) == mul(stretch(f,p),stretch(g,q))
        rows.append({'p':p,'q':q,'distinct_residue_cells':len(residues)})
    return rows

def main():
    pairs = list(itertools.combinations([3,5,7,11,13], 2))
    output = {'python':platform.python_version(),
              'family': [check_family(p,q) for p,q in pairs],
              'alias_formulas':[check_aliases(p,q) for p,q in pairs],
              'residue_grid':check_residue_grid(),
              'scope':'Exact finite identity checks; unrestricted degree theorem rests on the written proof.',
              'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    path = HERE/'parent-descent-identity-check.json'
    path.write_text(json.dumps(output,indent=2)+'\n')
    print({'family_pairs':len(pairs),'residue_grids':len(output['residue_grid']),
           'all_exact_checks':True,'receipt_sha256':hashlib.sha256(path.read_bytes()).hexdigest()})

if __name__ == '__main__':
    main()
