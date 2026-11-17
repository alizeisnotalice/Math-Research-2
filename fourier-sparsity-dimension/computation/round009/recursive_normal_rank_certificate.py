#!/usr/bin/env python3
"""Generate and replay exact finite certificates for recursive normal-rank bounds.

Generic certificates use exact normal-equation inconsistency, without point
enumeration.  The exhaustive n=4 audit additionally reconstructs all sixteen
points independently to cross-check that predicate and extension packing.
"""
import argparse, hashlib, itertools, json, random, sys, time
from fractions import Fraction

def rank(rows, n):
    a=[x & ((1<<n)-1) for x in rows if x & ((1<<n)-1)]; r=0
    for b in range(n-1,-1,-1):
        p=next((i for i in range(r,len(a)) if a[i]>>b&1),None)
        if p is None: continue
        a[r],a[p]=a[p],a[r]
        for i in range(len(a)):
            if i!=r and a[i]>>b&1:a[i]^=a[r]
        r+=1
    return r

def canonical(rows,n):
    a=[x & ((1<<n)-1) for x in rows if x & ((1<<n)-1)]; r=0
    for b in range(n-1,-1,-1):
        p=next((i for i in range(r,len(a)) if a[i]>>b&1),None)
        if p is None: continue
        a[r],a[p]=a[p],a[r]
        for i in range(len(a)):
            if i!=r and a[i]>>b&1:a[i]^=a[r]
        r+=1
    return tuple(a[:r])

def in_span(x, rows, n): return rank(list(rows)+[x],n)==rank(rows,n)
def span_values(rows):
    """Each vector in this tiny row span maps to its coefficient parity mask."""
    d=len(rows); out={}
    for q in range(1<<d):
        v=0
        for i in range(d):
            if q>>i&1:v^=rows[i]
        out[v]=q
    return out
def val_on(v, rows, rhs):
    q=span_values(rows)[v]
    return sum(((q>>i)&1)*rhs[i] for i in range(len(rows)))&1

def independent_system(rows, rhs, n):
    # Equation rows*x=rhs is consistent exactly when augmenting cannot raise rank.
    return rank(rows,n)==rank([a | (int(b)<<n) for a,b in zip(rows,rhs)],n+1)

def pair_disjoint(a,b,n):
    return not independent_system(a['rows']+b['rows'],a['rhs']+b['rhs'],n)

def all_disjoint(family,n):
    return all(pair_disjoint(family[i],family[j],n)
               for i in range(len(family)) for j in range(i))

def well_formed(f,n,c=None):
    return (len(f['rows'])==len(f['rhs']) and
            all(isinstance(x,int) and 0 < x < (1<<n) for x in f['rows']) and
            all(b in (0,1) for b in f['rhs']) and
            rank(f['rows'],n)==len(f['rows']) and
            (c is None or len(f['rows'])<=c))

def point_mask_small(f,n):
    """Independent brute-force predicate, only used for n=4 audit checks."""
    return sum(1<<x for x in range(1<<n)
               if all((bin(a&x).count('1')&1)==b for a,b in zip(f['rows'],f['rhs'])))

def intersection_restriction(f, other_rows, n):
    """Return an independent basis for W_f intersect span(other), with its rhs."""
    vals=[]
    for v in span_values(f['rows']):
        if v and in_span(v,other_rows,n): vals.append(v)
    basis=[]
    for v in vals:
        if rank(basis+[v],n)>len(basis):basis.append(v)
    return {'rows':basis,'rhs':[val_on(v,f['rows'],f['rhs']) for v in basis]}

def flat(rows,rhs,label=''):
    assert len(rows)==len(rhs)
    return {'rows':list(rows),'rhs':list(rhs),'label':label}

def encode_flat(f): return {'normal_rows':f['rows'],'rhs':f['rhs'],'label':f.get('label','')}
def decode_flat(f): return flat(f['normal_rows'],f['rhs'],f.get('label',''))

def minimal_rank_core(family,n):
    """Deletion certificate: retain an inclusion-minimal subfamily with same span."""
    keep=list(range(len(family))); target=rank([x for f in family for x in f['rows']],n); deleted=[]
    changed=True
    while changed:
        changed=False
        for pos,i in enumerate(keep[:]):
            trial=[j for j in keep if j!=i]
            rt=rank([x for j in trial for x in family[j]['rows']],n)
            if rt==target:
                keep.remove(i);deleted.append(i);changed=True;break
    assert rank([x for i in keep for x in family[i]['rows']],n)==target
    return keep,deleted,target

def build_certificate(family,n,c,depth=0):
    assert all(well_formed(f,n,c) for f in family) and all_disjoint(family,n)
    keep,deleted,total_rank=minimal_rank_core(family,n)
    core=[family[i] for i in keep]
    private=[]; extension=[]
    for i,f in enumerate(core):
        others=[x for j,g in enumerate(core) if j!=i for x in g['rows']]
        u=intersection_restriction(f,others,n)
        pi=len(f['rows'])-len(u['rows'])
        assert pi>=1
        private.append({'flat':encode_flat(f),'U':encode_flat(u),'p':pi})
        extension.append(flat(u['rows'],u['rhs'],f"U({f.get('label','')})"))
    assert all_disjoint(extension,n)
    kraft=sum(1<<x['p'] for x in private)
    assert kraft <= (1<<c)
    extension_mass=sum((Fraction(1,1<<len(f['rows'])) for f in extension),Fraction())
    assert extension_mass<=1
    sum_p=sum(x['p'] for x in private)
    assert sum_p<= (1<<(c-1)) if c else sum_p==0
    rank_u=rank([x for f in extension for x in f['rows']],n)
    assert total_rank <= sum_p+rank_u
    child=None
    if c:
        child=build_certificate(extension,n,c-1,depth+1) if extension else None
        child_bound=child['claimed_bound'] if child else 0
    else:
        assert rank_u==0; child_bound=0
    claimed=(1<<c)-1
    assert total_rank<=claimed
    return {'n':n,'c':c,'depth':depth,'input': [encode_flat(f) for f in family],
            'input_rank':rank([x for f in family for x in f['rows']],n),
            'target_rank':total_rank,'kept_input_indices':keep,'deleted_rank_redundant_indices':deleted,
            'core':[encode_flat(f) for f in core],'private':private,
            'extensions_pairwise_disjoint':True,'kraft_sum_2powp':kraft,'kraft_bound_2powc':1<<c,
            'extension_measure_sum':str(extension_mass),'extension_measure_bound':'1',
            'sum_private_dimensions':sum_p,'rank_sum_U':rank_u,
            'span_dimension_upper':sum_p+rank_u,'child_bound':child_bound,
            'claimed_bound':claimed,'child':child}

def verify_certificate(cert):
    n,c=cert['n'],cert['c']; fam=[decode_flat(f) for f in cert['input']]
    assert isinstance(n,int) and n>=1 and isinstance(c,int) and c>=0
    assert all(well_formed(f,n,c) for f in fam) and all_disjoint(fam,n)
    assert cert['input_rank']==rank([x for f in fam for x in f['rows']],n)
    assert cert['target_rank']==cert['input_rank']
    keep=cert['kept_input_indices']; deleted=cert['deleted_rank_redundant_indices']
    assert sorted(keep+deleted)==list(range(len(fam))) and len(set(keep))==len(keep)
    core=[fam[i] for i in keep]; target=cert['target_rank']
    assert cert['core']==[encode_flat(f) for f in core]
    assert len(cert['private'])==len(core)
    assert rank([x for f in core for x in f['rows']],n)==target
    # Each deleted index was redundant at deletion-time in the documented order is
    # stronger than needed only if we replay in forward order from original set.
    alive=list(range(len(fam)))
    for i in deleted:
        assert i in alive; trial=[j for j in alive if j!=i]
        assert rank([x for j in trial for x in fam[j]['rows']],n)==target; alive=trial
    assert alive==keep
    ps=[]; ext=[]
    for i,f in enumerate(core):
        others=[x for j,g in enumerate(core) if j!=i for x in g['rows']]
        u=intersection_restriction(f,others,n); got=cert['private'][i]
        assert encode_flat(f)==got['flat'] and encode_flat(u)==got['U']
        p=len(f['rows'])-len(u['rows']); assert p==got['p'] and p>=1
        ps.append(p); ext.append(flat(u['rows'],u['rhs'],f"U({f.get('label','')})"))
    assert all_disjoint(ext,n) and cert['extensions_pairwise_disjoint']
    kraft=sum(1<<p for p in ps); assert kraft==cert['kraft_sum_2powp']<=1<<c
    extension_mass=sum((Fraction(1,1<<len(f['rows'])) for f in ext),Fraction())
    assert str(extension_mass)==cert['extension_measure_sum'] and extension_mass<=1
    su=sum(ps); ru=rank([x for f in ext for x in f['rows']],n)
    assert su <= ((1<<(c-1)) if c else 0)
    assert su==cert['sum_private_dimensions'] and ru==cert['rank_sum_U']
    assert target<=su+ru==cert['span_dimension_upper']
    if c==0:
        assert cert['child'] is None and ru==0
    else:
        if not ext:
            assert cert['child'] is None and ru==0 and cert['child_bound']==0
        else:
            assert cert['child'] is not None; verify_certificate(cert['child'])
            assert cert['child']['n']==n and cert['child']['c']==c-1
            assert cert['child']['target_rank']==cert['child']['input_rank']
            assert cert['child']['input']==[encode_flat(f) for f in ext]
            assert ru<=cert['child']['claimed_bound']
    assert target<=cert['claimed_bound']==(1<<c)-1

def verify_point_certificate_n4(cert):
    """A separate truth-mask audit for the n=4 exhaustive cases only."""
    assert cert['n']==4
    fam=[decode_flat(f) for f in cert['input']]
    masks=[point_mask_small(f,4) for f in fam]
    for i in range(len(fam)):
        for j in range(i):
            assert pair_disjoint(fam[i],fam[j],4)==(not bool(masks[i]&masks[j]))
    ext=[decode_flat(x['U']) for x in cert['private']]
    emasks=[point_mask_small(f,4) for f in ext]
    assert all(not (emasks[i]&emasks[j]) for i in range(len(emasks)) for j in range(i))
    assert sum(bin(m).count('1') for m in emasks)<=16
    if cert['child'] is not None: verify_point_certificate_n4(cert['child'])

def all_subspaces(n,d):
    spaces={()}
    for _ in range(d):
        nxt=set()
        for s in spaces:
            for x in range(1,1<<n):
                if rank(list(s)+[x],n)>len(s):nxt.add(canonical(list(s)+[x],n))
        spaces=nxt
    return sorted(spaces)

def all_n4_flats():
    out=[]
    for d in range(3):
        for s in all_subspaces(4,d):
            for rhs in itertools.product(range(2),repeat=d):out.append(flat(s,list(rhs),f'c{d}:{s}:{rhs}'))
    return out

def exhaustive_n4_c2():
    flats=all_n4_flats(); assert len(flats)==171
    masks=[point_mask_small(f,4) for f in flats]
    compatible=[]
    for i in range(len(flats)):
        compatible.append({j for j in range(i+1,len(flats)) if pair_disjoint(flats[i],flats[j],4)})
        for j in range(i): assert pair_disjoint(flats[i],flats[j],4)==(not bool(masks[i]&masks[j]))
    counts={}; maxrank=0; sample=None; tested=0; rolling=bytes(32)
    def dfs(chosen,candidates):
        nonlocal maxrank,sample,tested,rolling
        cert=build_certificate([flats[i] for i in chosen],4,2); verify_certificate(cert); verify_point_certificate_n4(cert)
        payload=json.dumps(cert,sort_keys=True,separators=(',',':')).encode()
        rolling=hashlib.sha256(rolling+payload).digest()
        r=cert['target_rank']; tested+=1; counts[len(chosen)]=counts.get(len(chosen),0)+1
        if r>maxrank:maxrank=r;sample=cert
        for q,i in enumerate(sorted(candidates)):
            rest={j for j in candidates if j>i and j in compatible[i]}
            dfs(chosen+[i],rest)
    dfs([],set(range(len(flats))))
    return {'flat_count':len(flats),'families_including_empty':tested,'family_counts_by_size':counts,
            'max_normal_rank':maxrank,'case_certificate_rolling_sha256':rolling.hex(),
            'sample_maxrank_certificate':sample}

def random_independent_vectors(rng,n,count):
    out=[]
    while len(out)<count:
        x=rng.randrange(1,1<<n)
        if rank(out+[x],n)>len(out):out.append(x)
    return out

def mixed_tree(rng,n=31,maxdepth=5):
    queries=random_independent_vectors(rng,n,(1<<maxdepth)-1); leaves=[]; nodes=[]
    def rec(node,depth,rows,rhs):
        # no stopping before depth 2; then deterministic random mixed leaves.
        if depth==maxdepth or (depth>=2 and rng.random()<0.36):
            leaves.append(flat(rows,rhs,f'leaf{len(leaves)}'));return
        q=queries[node]; nodes.append({'node':node,'query':q,'depth':depth})
        rec(2*node+1,depth+1,rows+[q],rhs+[0]);rec(2*node+2,depth+1,rows+[q],rhs+[1])
    rec(0,0,[],[])
    return nodes,leaves

def tree_trials(seed=2026111709):
    rng=random.Random(seed); trials=[]
    for t in range(10):
        nodes,leaves=mixed_tree(rng)
        # fixed, reconstitutable pseudo-random subfamilies; all parent leaves are disjoint.
        for s in range(12):
            subset=[f for f in leaves if rng.getrandbits(1)]
            cert=build_certificate(subset,31,5); verify_certificate(cert)
            trials.append({'tree_index':t,'subset_index':s,'node_queries':nodes,
                           'leaf_count_total':len(leaves),'selected_leaf_indices':[leaves.index(f) for f in subset],
                           'certificate':cert})
    return {'seed':seed,'ambient_n':31,'max_depth':5,'trial_count':len(trials),'trials':trials}

def full_tree_controls():
    """Equality controls: complete depth-c trees, and trees with one leaf removed."""
    out=[]
    for c in range(1,6):
        q=[1<<i for i in range((1<<c)-1)]
        leaves=[]
        def rec(node,depth,rows,rhs):
            if depth==c:
                leaves.append(flat(rows,rhs,f'd{c}-leaf{len(leaves)}')); return
            rec(2*node+1,depth+1,rows+[q[node]],rhs+[0])
            rec(2*node+2,depth+1,rows+[q[node]],rhs+[1])
        rec(0,0,[],[])
        for name,fam in [('full',leaves),('delete_leaf_0',leaves[1:])]:
            cert=build_certificate(fam,31,c); verify_certificate(cert)
            assert cert['target_rank']==(1<<c)-1
            out.append({'c':c,'case':name,'query_normals':q,
                        'input':[encode_flat(f) for f in fam],'certificate':cert})
    return {'ambient_n':31,'controls':out}

def remote_packing_witnesses(path):
    raw=json.load(open(path)); result=raw['result']; records=[]
    for i,x in enumerate(result['pilot_records']): records.append((f'pilot-{i}',x['result']))
    for i,x in enumerate(result['jobs']): records.append((f'job-{i:02d}',x['result']))
    checked=[]
    for label,r in records:
        n=r['n']; fam=[flat(f['normal_rows'],f['rhs'],label) for f in r['flats']]
        cert=build_certificate(fam,n,3); verify_certificate(cert)
        assert cert['target_rank']==r['best_normal_rank'] and all_disjoint(fam,n)
        checked.append({'label':label,'n':n,'seed':r['seed'],'iterations':r['iterations_executed'],
                        'input':[encode_flat(f) for f in fam],'certificate':cert})
    return {'receipt_path':path,'receipt_sha256':hashlib.sha256(open(path,'rb').read()).hexdigest(),
            'search_source_sha256':result['source_sha256'],'records':checked}

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--out',required=True,help='output JSON path')
    ap.add_argument('--tree-seed',type=int,default=2026111709,help='fixed seed for mixed tree trials')
    ap.add_argument('--remote-packing',default='/Volumes/data/CodexTemp/T2/a28/r008/remote-packing-execution.json',help='r008 14-worker receipt JSON')
    a=ap.parse_args()
    began=time.perf_counter(); exact=exhaustive_n4_c2(); trees=tree_trials(a.tree_seed); controls=full_tree_controls(); remote=remote_packing_witnesses(a.remote_packing)
    source=hashlib.sha256(open(__file__,'rb').read()).hexdigest()
    out={'status':'passed','method':'exact_GF2_recursive_normal_rank_certificate',
         'n4_c_le_2_exhaustive':exact,'mixed_depth_tree_trials':trees,'full_tree_equality_controls':controls,
         'r008_remote_packing_witnesses':remote,
         'source_sha256':source,'python':sys.version,'elapsed_seconds':time.perf_counter()-began}
    open(a.out,'w').write(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'status':'passed','n4_families':exact['families_including_empty'],
                      'n4_max_rank':exact['max_normal_rank'],'tree_trials':trees['trial_count'],
                      'remote_records':len(remote['records']),
                      'elapsed_seconds':out['elapsed_seconds']}))
if __name__=='__main__':main()
