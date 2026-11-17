#!/usr/bin/env python3
"""Independent exact finite audit of the seven recorded r007 dim-15 extensions."""
import argparse, hashlib, json, platform, sys, time
from pathlib import Path

N, P, FULL = 6, 64, (1 << 64) - 1

def echelon(vectors):
    piv={}
    for x in vectors:
        while x:
            b=x.bit_length()-1
            if b in piv: x ^= piv[b]
            else: piv[b]=x; break
    return piv
def rank(vectors): return len(echelon(vectors))
def contains(x,piv):
    while x:
        b=x.bit_length()-1
        if b not in piv:return False
        x ^= piv[b]
    return True
def labeled_echelon(vectors):
    piv={}
    for j,x in enumerate(vectors):
        lab=1<<j
        while x:
            b=x.bit_length()-1
            if b in piv: x^=piv[b][0]; lab^=piv[b][1]
            else:piv[b]=(x,lab);break
    return piv
def coord(x,piv):
    lab=0
    while x:
        b=x.bit_length()-1
        if b not in piv:raise AssertionError('not in span')
        x^=piv[b][0];lab^=piv[b][1]
    return lab
def span(gens):
    out=[0]
    for g in gens: out += [x^g for x in out]
    return out
def fwht(a):
    a=a[:]; step=1
    while step<len(a):
        for s in range(0,len(a),2*step):
            for j in range(s,s+step): a[j],a[j+step]=a[j]+a[j+step],a[j]-a[j+step]
        step*=2
    return a
def q_mask(permutation,h_bits):
    return sum((((((x&7)&permutation[x>>3]).bit_count() ^ ((h_bits>>(x>>3))&1))&1)<<x) for x in range(P))
def bent(mask): return all(abs(x)==8 for x in fwht([-1 if mask>>x&1 else 1 for x in range(P)]))
def affine_block(values,block):
    pts=[x for x in range(P) if block>>x&1]; base=pts[0]; dirs=[x^base for x in pts]
    if len(pts)!=8 or rank(dirs)!=3:return False,[]
    basis=list(echelon(dirs).values())
    for d in basis:
        slope=values[base^d]^values[base]
        if any(values[x^d]^values[x]!=slope for x in pts):return False,basis
    return True,basis
def c_space(): return [FULL]+[sum(((x>>j)&1)<<x for x in range(P)) for j in range(N)]
def minword_structure(words):
    minimum=min(x.bit_count() for x in words if x); ws=[x for x in words if x.bit_count()==minimum]
    if minimum != 8: raise AssertionError('minimum weight')
    adj={i:{j for j,y in enumerate(ws) if i!=j and (ws[i]&y).bit_count()==4} for i in range(len(ws))}
    if any((ws[i]&ws[j]).bit_count() not in (0,4) for i in range(len(ws)) for j in range(i)):raise AssertionError('overlap')
    seen=set(); comps=[]
    for i in range(len(ws)):
        if i not in seen:
            todo=[i];seen.add(i); comp=[]
            while todo:
                u=todo.pop();comp.append(u)
                for v in adj[u]:
                    if v not in seen:seen.add(v);todo.append(v)
            comps.append(comp)
    six=[c for c in comps if len(c)==6]; isolated=[c for c in comps if len(c)==1]
    if len(six)+len(isolated)!=len(comps):raise AssertionError('component size')
    if any(any(len(adj[i]&set(c))!=4 for i in c) for c in six):raise AssertionError('E4 degree')
    index={x:i for i,x in enumerate(ws)}
    for i in range(len(ws)):
        for j in adj[i]:
            if ws[i]^ws[j] not in index:raise AssertionError('triangle A4 closure')
    used = 0
    for comp in comps:
        union = 0
        for i in comp: union |= ws[i]
        if used & union: raise AssertionError('component supports overlap')
        used |= union
        if len(comp) == 1:
            if union.bit_count() != 8: raise AssertionError('isolated size')
            continue
        atoms = {}
        for x in range(P):
            label = sum(((ws[i] >> x) & 1) << j for j,i in enumerate(comp))
            if label: atoms[label] = atoms.get(label, 0) | (1 << x)
        atoms = list(atoms.values())
        if len(atoms) != 4 or any(a.bit_count() != 4 for a in atoms):
            raise AssertionError('four half-flat atoms')
        directions, bases = [], []
        for atom in atoms:
            pts = [x for x in range(P) if atom >> x & 1]
            base = pts[0]; dirs = {x ^ base for x in pts}
            if rank(dirs) != 2: raise AssertionError('atom affine dimension')
            directions.append(dirs); bases.append(base)
        if any(d != directions[0] for d in directions): raise AssertionError('atom direction')
        if bases[0] ^ bases[1] ^ bases[2] ^ bases[3] not in directions[0]:
            raise AssertionError('fixed fourth coset')
        edges = {a ^ b for j,a in enumerate(atoms) for b in atoms[:j]}
        if edges != {ws[i] for i in comp} or rank(edges) != 3:
            raise AssertionError('E4 equality')
    return ws,adj,six,isolated

ap=argparse.ArgumentParser()
ap.add_argument('--candidates',type=Path,required=True);ap.add_argument('--extensions',type=Path,required=True);ap.add_argument('--audit',type=Path,required=True);ap.add_argument('--out',type=Path,required=True)
a=ap.parse_args(); began=time.monotonic()
rc,re,ra=a.candidates.read_bytes(),a.extensions.read_bytes(),a.audit.read_bytes()
CAND,EXT,AUD=json.loads(rc),json.loads(re),json.loads(ra)
seed=next(r for r in AUD['rows'] if r['source_file'].endswith('corrected-results.json') and r['row_index']==1)
q=q_mask(seed['permutation'],seed['h_bits']); W=CAND['W_generators']; K=set(CAND['K'])
if rank(W)!=11 or not all(contains(c,echelon(W)) for c in c_space()):raise AssertionError('W')
if EXT['dimension15_extensions']!=7 or len(EXT['records'])!=7:raise AssertionError('record count')
rows=[]; partial=None
for n,rec in enumerate(EXT['records']):
    QB=rec['quotient_basis']; G=W+QB
    if len(QB)!=4 or rank(G)!=15 or not all(x in K for x in QB):raise AssertionError(('basis',n))
    blocks=[int(x,16) for x in rec['partition_masks_hex']]
    lp=labeled_echelon(G); labels=[coord(b,lp) for b in blocks]
    if any(b.bit_count()!=8 for b in blocks) or len(set(blocks))!=8:raise AssertionError('blocks')
    union=0
    for b in blocks:
        if union&b:raise AssertionError('not disjoint')
        union|=b
    if union!=FULL:raise AssertionError('not cover')
    directions=[]
    for b in blocks:
        funcs=[[q>>x&1 for x in range(P)]]+[[(g>>x)&1 for x in range(P)] for g in G]
        ok,dirs=affine_block(funcs[0],b)
        if not ok or any(not affine_block(f,b)[0] for f in funcs[1:]):raise AssertionError(('affine',n))
        directions.append(dirs)
    dual=[]
    for b in blocks:
        out=fwht([(-1 if q>>x&1 else 1) if b>>x&1 else 0 for x in range(P)])
        if any(x not in (-8,0,8) for x in out) or sum(x!=0 for x in out)!=8:raise AssertionError('dual size')
        dual.append(out)
    if any(sum(o[x]!=0 for o in dual)!=1 for x in range(P)):raise AssertionError('dual partition')
    normals=[]
    for b,dirs in zip(blocks,directions):
        base=(b&-b).bit_length()-1
        for d in dirs:
            slopes=sum((((g>>(base^d))^(g>>base))&1)<<j for j,g in enumerate(G))
            normals.append(slopes|(d<<15))
    jr, sr=rank(normals),rank([x&((1<<15)-1) for x in normals])
    if jr!=7 or sr!=7:raise AssertionError('normal rank')
    words=span(G); ws,adj,six,iso=minword_structure(words)
    h,p=len(six),len(iso)
    if rank(ws)!=3*h+p or 3*h+p>12:raise AssertionError('min rank')
    if not six:raise AssertionError('no E4')
    if not all(bent(q^x) for x in words): raise AssertionError('non-bent full slice')
    rows.append({'record_index':n,'rank_G':rank(G),'all_32768_slices_bent':True,
                 'partition_block_coordinates_hex':[hex(x) for x in labels],'affine_blocks':8,'dual_blocks_partition':True,
                 'joint_normal_rank':jr,'slope_rank':sr,'minimum_weight':8,'minimum_words':len(ws),
                 'E4_components':h,'isolated_minimum_words':p,'minimum_span_rank':rank(ws),'triangle_A4_fixed':True,'E4_contained':True})
    if partial is None:
        comp=six[0]; pair=next((i,j) for i in comp for j in adj[i] if j in comp and i<j)
        x,y=ws[pair[0]],ws[pair[1]]
        z=next(ws[k] for k in comp if rank(c_space()+[x,y,ws[k]])==10)
        if rank(c_space()+[x,y])!=9 or rank(c_space()+[x,y,z])!=10:raise AssertionError('partial rank')
        partial={'record_index':n,'intersecting_minwords_hex':[hex(x),hex(y)],'completion_hex':hex(z),
                 'partial_dim':9,'completed_dim':10,'completion_raises_dim_by':1,
                 'all_8_completed_quotient_slices_bent':all(bent(q^t) for t in span([x,y,z]))}
        if not partial['all_8_completed_quotient_slices_bent']:raise AssertionError('partial slices')
result={'status':'passed','scope':'finite verification of supplied r007 seven extensions; no general proof claim',
        'candidate_sha256':hashlib.sha256(rc).hexdigest(),'extensions_sha256':hashlib.sha256(re).hexdigest(),'audit_sha256':hashlib.sha256(ra).hexdigest(),
        'seed_row':'corrected-results row 1','W_rank':11,'extension_count':7,'rows':rows,'partial_V_certificate':partial,
        'elapsed_seconds':time.monotonic()-began,'python':sys.version,'platform':platform.platform(),'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
a.out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({'status':'passed','extensions':7,'partial_completed_slices':8,'elapsed_seconds':result['elapsed_seconds']}))
