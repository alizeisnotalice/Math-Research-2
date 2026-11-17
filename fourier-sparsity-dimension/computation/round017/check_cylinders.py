#!/usr/bin/env python3
"""Independent exact audit of round-017 joint parameter--frequency cylinders.

No project helper is imported.  Arithmetic is GF(2) bit arithmetic and
integer Walsh sums only.  It reconstructs local ANFs, freezes non-star
quadratic branches, constructs (7)--(8), and checks their geometry.
"""
from __future__ import annotations
import argparse, hashlib, json, platform, time
from pathlib import Path

HERE=Path(__file__).resolve().parent
DEFAULT_INPUT=HERE/'nonparallel-cases.json'
DEFAULT_OUTPUT=HERE/'cylinder-audit.json'

def dot(a,b): return (a&b).bit_count()&1
def rank(xs):
    piv={}
    for x in xs:
        while x:
            p=x.bit_length()-1
            if p in piv:x^=piv[p]
            else:piv[p]=x;break
    return len(piv)
def indep_basis(xs):
    out=[]
    for x in xs:
        if rank(out+[x])>len(out):out.append(x)
    return out
def span(bs):
    out={0}
    for b in bs:out|={x^b for x in tuple(out)}
    return out
def xor_select(xs,mask):
    out=0
    for j,x in enumerate(xs):
        if mask>>j&1:out^=x
    return out
def nullspace(rows,n):
    return [x for x in range(1<<n) if all(not dot(r,x) for r in rows)]
def basis_of_set(S): return indep_basis(sorted(S,reverse=True))
def affine_info(points,n):
    assert points
    p=points[0]; db=indep_basis([x^p for x in points[1:]])
    full={p^x for x in span(db)}
    return {'size':len(points),'direction_dimension':len(db),'is_affine':full==set(points),'base':p,'direction_basis':db},db
def orthogonal_basis(ds,n):
    return basis_of_set(nullspace(ds,n))
def bit(mask,x):return (mask>>x)&1
def phase(mask,pts):return [bit(mask,x) for x in pts]
def mobius(vals,d):
    a=vals[:]
    for j in range(d):
        for m in range(1<<d):
            if m>>j&1:a[m]^=a[m^(1<<j)]
    return a
def local_points(p,hb):
    return [p^__import__('functools').reduce(lambda a,b:a^b,
            (hb[j] for j in range(len(hb)) if m>>j&1),0) for m in range(1<<len(hb))]
def local_anf(mask,pts):return mobius(phase(mask,pts),(len(pts)).bit_length()-1)
def qvec(anf,d):
    out=[]
    for i in range(d):
        for j in range(i+1,d):out.append(anf[(1<<i)|(1<<j)])
    return sum(v<<i for i,v in enumerate(out))
def linvec(anf,d):return sum(anf[1<<i]<<i for i in range(d))
def prod_qvec(a,b,d):
    vals=[]
    for i in range(d):
        for j in range(i+1,d): vals.append(((a>>i&1)&(b>>j&1))^((a>>j&1)&(b>>i&1)))
    return sum(v<<i for i,v in enumerate(vals))
def choose_common_factor(qs,d):
    # accepts zero only if never used here; all double blocks rank two.
    good=[]
    for a in range(1,1<<d):
        if all(any(prod_qvec(a,b,d)==q for b in range(1<<d)) for q in qs):good.append(a)
    return min(good) if good else None
def choose_b(a,q,d):
    pivot=(a&-a).bit_length()-1
    bs=[b for b in range(1<<d) if not (b>>pivot&1) and prod_qvec(a,b,d)==q]
    assert len(bs)==1,(a,q,bs)
    return bs[0]
def walsh_local(mask,pts):
    d=(len(pts)).bit_length()-1
    return [sum((1 if bit(mask,x)==0 else -1)* (1 if dot(l,m)==0 else -1)
                for m,x in enumerate(pts)) for l in range(1<<d)]
def restrict_freq(xi,hb):return sum(dot(xi,h)<<j for j,h in enumerate(hb))
def partial_walsh(mask,block,xi):return sum((1 if bit(mask,x)==0 else -1)*(1 if dot(xi,x)==0 else -1) for x in block)

def analyze(case):
    d=len(case['basis_masks']); bm=case['basis_masks']; f0=case['f0_mask']; blocks=case['blocks']
    assert f0>>64==0 and all(m>>64==0 for m in bm)
    assert rank(bm)==d
    flat=[x for A in blocks for x in A]
    assert sorted(flat)==list(range(64)) and len(flat)==len(set(flat))
    desc=[]
    for A in blocks:
        p=A[0]; H={x^p for x in A}; hb=basis_of_set(H)
        assert span(hb)==H
        pts=local_points(p,hb); assert set(pts)==set(A)
        desc.append((p,hb,pts))
    double=[i for i,A in enumerate(blocks) if len(A)==16]
    minimal=[i for i,A in enumerate(blocks) if len(A)==8]
    # Local quadratic columns for full parameter code.
    q0={}; qcols={}; star={}; common={}
    for i in double:
        _,hb,pts=desc[i]; ld=len(hb)
        a0=local_anf(f0,pts); q0[i]=qvec(a0,ld)
        cols=[qvec(local_anf(f0^m,pts),ld)^q0[i] for m in bm]
        qcols[i]=cols
        qs=[q0[i]^xor_select(cols,p) for p in range(1<<d)]
        aa=choose_common_factor(qs,ld)
        star[i]=aa is not None; common[i]=aa
        assert all(q!=0 for q in qs)
    nonstar=[i for i in double if not star[i]]
    # freezing rows (coordinate functionals in parameter space)
    rows=[]
    for i in nonstar:
        outbits=max(1,max(qcols[i]).bit_length())
        for c in range(outbits):
            row=0
            for j in range(d):
                if qcols[i][j]>>c&1:row|=1<<j
            rows.append(row)
    rows=[r for r in rows if r]
    rfreeze=rank(rows); assert rfreeze<=3*len(nonstar)
    v0all=nullspace(rows,d); v0basis=indep_basis(v0all); d0=len(v0basis)
    assert d0==d-rfreeze
    assert len(v0all)==1<<d0
    def op(t):
        x=0
        for k in range(d0):
            if t>>k&1:x^=v0basis[k]
        return x
    def ph(t):
        p=op(t)
        x=f0
        for j in range(d):
            if p>>j&1:x^=bm[j]
        return x
    # Recover a,b,ell forms on frozen code.
    forms={}
    for i in double:
        _,hb,pts=desc[i]; ld=len(hb)
        qs=[qvec(local_anf(ph(t),pts),ld) for t in range(1<<d0)]
        aa=choose_common_factor(qs,ld)
        if aa is None:
            raise AssertionError({'case':case['name'],'block':i,'V0_dim':d0,'quadratics':sorted(set(qs)),
                                  'full_star':star[i],'full_common':common[i],'freeze_rows':rows})
        dat=[]
        for t in range(1<<d0):
            an=local_anf(ph(t),pts); q=qvec(an,ld); bb=choose_b(aa,q,ld)
            ell=linvec(an,ld)^(aa&bb)
            # exact polynomial reconstruction check (constant ignored)
            for m in range(1<<ld):
                val=(an[0]^dot(ell,m)^ (dot(aa,m)&dot(bb,m)))
                assert val==phase(ph(t),pts)[m]
            assert bb not in (0,aa)
            dat.append((ell,bb))
        forms[i]=(aa,dat)
    for i in minimal:
        _,hb,pts=desc[i]; ld=len(hb); assert ld==3
        dat=[]
        for t in range(1<<d0):
            an=local_anf(ph(t),pts)
            assert all(not an[m] or m.bit_count()<=1 for m in range(1<<ld))
            assert qvec(an,ld)==0
            dat.append(linvec(an,ld))
        forms[i]=dat
    # Construct cylinders in parameter x frequency ambient.
    N=d0+6; cylinders=[]; labels=[]
    for i in double:
        _,hb,_=desc[i]; aa,dat=forms[i]
        for u in (0,1):
            pts=[]
            for t in range(1<<d0):
                ell,bb=dat[t]
                for xi in range(64):
                    rr=restrict_freq(xi,hb)^ell^(bb if u else 0)
                    if rr in (0,aa):pts.append(t|(xi<<d0))
            info,db=affine_info(pts,N); assert info['is_affine'] and info['direction_dimension']==N-3
            cylinders.append((pts,db));labels.append(f'double_{i}_u{u}')
    for i in minimal:
        _,hb,_=desc[i]; dat=forms[i];pts=[]
        for t in range(1<<d0):
            for xi in range(64):
                if restrict_freq(xi,hb)==dat[t]:pts.append(t|(xi<<d0))
        info,db=affine_info(pts,N);assert info['is_affine'] and info['direction_dimension']==N-3
        cylinders.append((pts,db));labels.append(f'minimal_{i}')
    sets=[set(a) for a,_ in cylinders]; union=set()
    for s in sets:
        assert not union&s;union|=s
    assert union==set(range(1<<N))
    normals=[]
    for _,db in cylinders:normals+=orthogonal_basis(db,N)
    nrank=rank(normals); assert nrank<=7
    proj=[x&((1<<d0)-1) for x in normals]; prank=rank(proj)
    kvals=nullspace(proj,d0); kbasis=indep_basis(kvals)
    assert len(kbasis)+prank==d0
    # Exact cylinder slices versus actual block Walsh support.
    for t in range(1<<d0):
        mask=ph(t)
        for bi,A in enumerate(blocks):
            actual={xi for xi in range(64) if partial_walsh(mask,A,xi)!=0}
            assert all(abs(partial_walsh(mask,A,xi))==8 for xi in actual)
            expected=set()
            for ci,lab in enumerate(labels):
                if (lab.startswith(f'double_{bi}_') or lab==f'minimal_{bi}'):
                    expected|={p>>d0 for p in sets[ci] if (p&((1<<d0)-1))==t}
            assert actual==expected,(case['name'],t,bi,len(actual),len(expected))
    # Restriction ranks for K directions (not including f0).
    kranks=[]
    for bi,A in enumerate(blocks):
        vs=[]
        for t in kbasis:
            p=op(t); vm=0
            for j in range(d):
                if p>>j&1:vm^=bm[j]
            vs.append(sum(bit(vm,x)<<q for q,x in enumerate(A)))
        kranks.append(rank(vs))
    concat=[]
    for t in kbasis:
        p=op(t); vm=0
        for j in range(d):
            if p>>j&1:vm^=bm[j]
        off=0;cc=0
        for A in blocks:
            cc|=sum(bit(vm,x)<<(off+q) for q,x in enumerate(A));off+=len(A)
        concat.append(cc)
    image_rank=rank(concat); local_limits=[2 if len(A)==16 else 1 for A in blocks]
    assert all(a<=b for a,b in zip(kranks,local_limits))
    assert image_rank<=sum(local_limits)
    return {
      'name':case['name'],'V_dimension':d,'freeze_nonstar_blocks':nonstar,
      'star_blocks':[i for i in double if star[i]],'freeze_rank_r':rfreeze,'V0_dimension':d0,
      'V0_phase_count':1<<d0,'common_factors_intrinsic':{str(i):forms[i][0] for i in double},
      'cylinder_count':len(cylinders),'cylinder_labels':labels,'joint_ambient_dimension':N,
      'cylinder_codimension':3,'disjoint_and_cover':True,
      'normal_union_rank':nrank,'parameter_projection_rank':prank,
      'K_dimension':len(kbasis),'kernel_restriction_ranks_by_block':kranks,
      'kernel_restriction_image_rank':image_rank,
      'kernel_bound_m1_plus_2m2':len(minimal)+2*len(double),
      'walsh_support_slices_exact':True,
    }

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--input',type=Path,default=DEFAULT_INPUT,
                    help='certificate JSON (default: script directory/nonparallel-cases.json)')
    ap.add_argument('--output',type=Path,default=DEFAULT_OUTPUT,
                    help='audit JSON (default: script directory/cylinder-audit.json)')
    ns=ap.parse_args(); start=time.monotonic()
    srcpath=ns.input.resolve(); outpath=ns.output.resolve()
    src=json.loads(srcpath.read_text()); out=[]
    for c in src['cases']:out.append(analyze(c))
    script=Path(__file__).resolve()
    result={'schema':'a28-r017-independent-cylinder-audit-v2','input':str(srcpath),
      'output':str(outpath),'input_sha256':hashlib.sha256(srcpath.read_bytes()).hexdigest(),
      'script_sha256':hashlib.sha256(script.read_bytes()).hexdigest(),'python':platform.python_version(),
      'arithmetic':'exact GF(2) and integer Walsh','elapsed_seconds':time.monotonic()-start,
      'cases':out,'failed_assertions':0,'exit_status':0}
    outpath.parent.mkdir(parents=True,exist_ok=True)
    outpath.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
