#!/usr/bin/env python3
"""Local exact audit of R19 non-star quotient-freeze cylinders.

This is deliberately only a four-variable alternating-form / cylinder model.
It does not construct a global common-bent code.  Its role is to check that
all h=0,1,2,3 entries in R19 §§4--5 are feasible and that the cylinder
equations, codimensions, and local K restriction budget agree exactly.
"""
import argparse, hashlib, itertools, json, multiprocessing as mp, os, platform, sys, time
from pathlib import Path

PAIRS=((0,1),(0,2),(0,3),(1,2),(1,3),(2,3))
PAIR_INDEX={p:i for i,p in enumerate(PAIRS)}

def dot(a,b): return (a & b).bit_count() & 1
def rank(vs):
    piv={}
    for x in vs:
        while x:
            p=x.bit_length()-1
            if p in piv: x^=piv[p]
            else: piv[p]=x; break
    return len(piv)
def span(bs):
    ans={0}
    for b in bs: ans|={x^b for x in tuple(ans)}
    return ans
def basis_of(S):
    out=[]; before=0
    for x in sorted(S, reverse=True):
        after=rank(out+[x])
        if after>before: out.append(x); before=after
    return out
def intersection_subspaces(sets):
    return set.intersection(*[set(x) for x in sets])

def col(form,j):
    ans=0
    for i in range(4):
        if i==j: continue
        p=(i,j) if i<j else (j,i)
        if form>>PAIR_INDEX[p]&1: ans|=1<<i
    return ans
def support_plane(form):
    return span([col(form,j) for j in range(4)])
def form_rank(form): return rank([col(form,j) for j in range(4)])
def exact_rank2(form): return form != 0 and form_rank(form)==2

def all_dim3_subspaces():
    found=set()
    for raw in range(1<<18):
        b=[(raw>>(6*j))&63 for j in range(3)]
        if rank(b)==3: found.add(tuple(sorted(span(b))))
    assert len(found)==1395
    return [set(x) for x in found]

def find_flag():
    """Find A0+L exact-rank2, non-star, and H1<H2<L with j=(1,1,0)."""
    r2=[x for x in range(64) if exact_rank2(x)]
    assert len(r2)==35
    for L in all_dim3_subspaces():
        if len(L)!=8: continue
        for a0 in r2:
            fam=[a0^x for x in L]
            if not all(exact_rank2(x) for x in fam): continue
            if intersection_subspaces([support_plane(x) for x in fam]) != {0}: continue
            nonzero=sorted(L-{0})
            for e1 in nonzero:
                h1=span([e1])
                if rank(intersection_subspaces([support_plane(a0^x) for x in h1]))!=1: continue
                for e2 in nonzero:
                    if rank([e1,e2])!=2: continue
                    h2=span([e1,e2])
                    if rank(intersection_subspaces([support_plane(a0^x) for x in h2]))!=1: continue
                    for e3 in nonzero:
                        if rank([e1,e2,e3])!=3: continue
                        assert rank(intersection_subspaces([support_plane(a0^x) for x in L]))==0
                        return a0,[e1,e2,e3]
    raise AssertionError('no required non-star flag found')

def encode(v,eta,beta,pdim,ny,h): return v | (eta<<pdim) | (beta<<(pdim+ny))
def affine_flat(points):
    assert points
    base=next(iter(points)); directions={x^base for x in points}
    r=rank(directions)
    return len(points)==(1<<r) and span(basis_of(directions))==directions, r
def nullspace(rows,n):
    return [x for x in range(1<<n) if all(dot(row,x)==0 for row in rows)]
def annihilator_of_direction(directions,n):
    return basis_of(nullspace(directions,n))

def audit_h(payload):
    h,a0,flag=payload
    H=span(flag[:h])
    E={gamma:support_plane(a0^gamma) for gamma in H}
    J=intersection_subspaces(list(E.values()))
    j=rank(J)
    expected_j=[2,1,1,0][h]
    assert j==expected_j
    # h quadratic-lift variables u, j affine-linear kernel variables w, and c.
    # For h=3 use the minimal actual Y dimension 5 and Z dimension 3; otherwise
    # use Y=F2^4,Z=F2^2.  Forms live in the first four Y coordinates.
    ny=5 if h==3 else 4
    pdim=h+j+1
    ambient=pdim+ny+h
    jbasis=basis_of(J)
    def B(coords):
        value=0
        for s in range(h):
            if coords>>s&1: value^=flag[s]
        return value
    cylinders=[]
    allnorm=[]
    for gamma in sorted(H):
        pts=set()
        for v in range(1<<pdim):
            u=v & ((1<<h)-1)
            w=(v>>h)&((1<<j)-1)
            ell=0
            for s,b in enumerate(jbasis):
                if w>>s&1: ell^=b
            for beta in range(1<<h):
                if B(beta^u) != gamma: continue # B(beta+T(v))=gamma
                for eta in range(1<<ny):
                    if (eta^ell) in E[gamma]:
                        pts.add(encode(v,eta,beta,pdim,ny,h))
        ok,d=affine_flat(pts)
        if not (ok and d==ambient-(h+(ny-2))):
            raise AssertionError({'h':h,'gamma':gamma,'flat':ok,'actual_dim':d,
                                  'expected_dim':ambient-(h+(ny-2)),'size':len(pts),
                                  'ambient':ambient})
        normals=annihilator_of_direction({x^next(iter(pts)) for x in pts},ambient)
        assert len(normals)==h+(ny-2)
        allnorm.extend(normals)
        cylinders.append({'gamma':gamma,'size':len(pts),'affine_direction_dimension':d,
                          'codimension':ambient-d,'normal_count':len(normals)})
    sets=[]
    for gamma in sorted(H):
        # recompute a compact point set to assert pairwise disjointness independently
        S=set()
        for v in range(1<<pdim):
            u=v&((1<<h)-1); w=(v>>h)&((1<<j)-1); ell=0
            for s,b in enumerate(jbasis):
                if w>>s&1: ell^=b
            for beta in range(1<<h):
                if B(beta^u)==gamma:
                    S|={encode(v,eta,beta,pdim,ny,h) for eta in range(1<<ny) if (eta^ell) in E[gamma]}
        sets.append(S)
    assert all(not(sets[i]&sets[z]) for i in range(len(sets)) for z in range(i))
    # Project normal vectors to parameter dual, then compute K=P^perp.
    pnorm=[n&((1<<pdim)-1) for n in allnorm]
    K=nullspace(pnorm,pdim)
    expected_K={v for v in range(1<<pdim) if (v&((1<<h)-1))==0}
    assert set(K)==expected_K
    phase_vectors=[]
    for v in K:
        w=(v>>h)&((1<<j)-1); c=(v>>(h+j))&1; ell=0
        for s,b in enumerate(jbasis):
            if w>>s&1: ell^=b
        phase_vectors.append(ell | (c<<ny)) # coefficients of affine phase on Y
    assert rank(phase_vectors)==1+j
    return {'h':h,'worker_pid':os.getpid(),'quadratic_frequency_dimension':h,
            'freeze_codimension':3-h,'j':j,'Y_dimension':ny,'Z_dual_dimension':h,
            'parameter_dimension':pdim,'ambient_dimension':ambient,'cylinder_count':len(cylinders),
            'cylinders':cylinders,'pairwise_disjoint':True,'parameter_projection_rank':rank(pnorm),
            'K_dimension':rank(K),'expected_K_dimension':1+j,
            'local_affine_phase_image_rank':rank(phase_vectors),'local_budget_bound':1+j}

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--output',type=Path,default=Path(__file__).with_name('nonstar-quotient-freeze-local-audit.json'))
    ap.add_argument('--workers',type=int,default=4)
    args=ap.parse_args(); assert 1<=args.workers<=4
    begun=time.monotonic(); a0,flag=find_flag()
    assert all(exact_rank2(a0^x) for x in span(flag))
    supports=[support_plane(a0^x) for x in span(flag)]
    assert intersection_subspaces(supports)=={0}
    with mp.Pool(processes=args.workers) as pool:
        rows=pool.map(audit_h,[(h,a0,flag) for h in range(4)])
    rows.sort(key=lambda x:x['h'])
    result={'schema':'a28-r019-nonstar-quotient-freeze-local-v1',
      'scope':('four-variable exact-rank2 alternating-form and quotient-freeze cylinder model only; '
               'not a global bent-code construction or proof of R19'),
      'arithmetic':'exact F2 enumeration','a0_form_hex':hex(a0),'joint_direction_basis_hex':[hex(x) for x in flag],
      'joint_family_size':8,'joint_family_all_exact_rank2':True,'joint_family_nonstar_common_support_dimension':0,
      'workers_requested':args.workers,'worker_pids':sorted({x['worker_pid'] for x in rows}),
      'rows':rows,'python':sys.version,'platform':platform.platform(),'elapsed_seconds':time.monotonic()-begun,
      'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'status':'passed','family':result['a0_form_hex'],
                      'basis':result['joint_direction_basis_hex'],'workers':result['worker_pids'],
                      'rows':[{k:r[k] for k in ['h','j','cylinder_count','K_dimension','local_affine_phase_image_rank']} for r in rows]},indent=2))
if __name__=='__main__': main()
