"""Certified Hausdorff radii covering the 501 five-body localization leaves."""
from argparse import ArgumentParser
from fractions import Fraction as F
from hashlib import sha256
from math import isqrt
from pathlib import Path
from io import BytesIO
import json
import sys

import numpy as np
from flint import arb, ctx

RAW = 2**44
UNIT = 2**28
AREA_DEN = 2**57
GLOBAL_LOW = F(120666205824033364, AREA_DEN)
PI_UP = F(22, 7)
BODY = ((3, -1), (5, 2), (7, 5), (9, 8))


def digest(payload):
    return sha256(payload).hexdigest()


def exact(x):
    x = F(x)
    return arb(x.numerator) / x.denominator


def checked(path, pin=None):
    payload = path.read_bytes()
    sha = digest(payload)
    if pin:
        assert sha == pin
    return json.loads(payload), sha


def ceil(x):
    return -(-x.numerator // x.denominator)


def sqrt_ceil(x):
    n = ceil(x)
    v = isqrt(n)
    return v + (v*v < n)


def enclose(x):
    q = round(float(x.mid()) * UNIT)
    assert arb(q-1) < UNIT*x < arb(q+1)
    return q-1, q+1


def read_leaves(main, tree):
    sys.path.insert(0, str(main))
    from tomar_leaf_box import leaf_box
    from verify_tomar_fan_table import PINNED_MANIFEST

    raw = (tree/'SHA256SUMS').read_bytes()
    manifest = digest(raw)
    assert manifest == PINNED_MANIFEST
    pins = {name.removeprefix('./'): h for h, name in (s.split() for s in raw.decode().splitlines())}
    high, hsha = checked(main/'five-global-localization.json',
                         'd85941b19b47cda83fc86268daa50569db2d80118b7fc722d4bf90adfa457504')
    clusters, csha = checked(main/'five-orientation-branches.json',
                             '654c72a37c28239680df4880a96319e23ed817ff9cc039b58f644bfc27935f2f')
    assert clusters['parent_sha256'] == hsha and clusters['input_manifest_sha256'] == manifest
    comps = clusters['components']
    assert sorted(map(len, comps)) == [25, 130, 346]
    comp = {tuple(key): group for group, keys in enumerate(comps) for key in keys}
    ids, low, up, lab = [], [], [], []
    for key, positions in sorted(high['retained'].items(), key=lambda row: int(row[0])):
        p = int(key)
        stem = f'part{p:04}'
        meta, _ = checked(tree/(stem+'.json'), pins[stem+'.json'])
        payload = (tree/(stem+'.npy')).read_bytes()
        assert digest(payload) == pins[stem+'.npy']
        codes = np.load(BytesIO(payload), allow_pickle=False)
        prefix = np.r_[1, 1+np.cumsum(np.where(codes < 11, 1, -1))]
        for pos in positions:
            assert 0 <= pos < len(codes) and codes[pos] >= 253
            c, w = leaf_box(prefix, codes, p, pos, meta)
            l, h = ([F(z)*RAW for z in arr] for arr in (c-w, c+w))
            assert all(v.denominator == 1 for v in l+h)
            ids.append([p, pos])
            low.append(list(map(int, l)))
            up.append(list(map(int, h)))
            lab.append(comp[(p, pos)])
    assert len(comp) == len(ids) == len({tuple(key) for key in ids}) == 501
    return (ids, np.array(low, np.int64), np.array(up, np.int64), np.array(lab),
            {'tree':manifest, 'outer_localization':hsha, 'orientation_branches':csha})


class Samples:
    def __init__(self, m):
        assert m >= 48 and m % 6 == 0
        self.m = m
        self.pi = arb.pi()
        self.angles = [2*self.pi*j/m for j in range(m)]
        self.radius = {n: 1/(2*(self.pi/(2*n)).cos()) for n, _ in BODY}
        self.lipschitz = {n:self.radius[n]*(self.pi/(2*n)).sin() for n, _ in BODY}
        self.normal = np.array([[enclose(a.cos()), enclose(a.sin())] for a in self.angles], np.int64)
        self.lattice = {}
        self.index = np.arange(m)
        self.perm = np.array([(sign*self.index+shift*m//3) % m
                              for sign in (1, -1) for shift in range(3)])

    def support(self, n, rho):
        rho = F(rho)
        key = (n, rho)
        if key in self.lattice:
            return self.lattice[key]
        gamma = self.pi/n
        r = self.radius[n]
        upper_tie = gamma+arb(1)/2**120
        base = exact(rho)
        period = float(2*gamma)
        values = []
        for i, u in enumerate(self.angles):
            raw = u-base
            choice = round((2*np.pi*i/self.m-float(rho))/period)
            folded = abs(raw-2*choice*gamma)
            assert folded < upper_tie
            if folded < gamma/2:
                h = r*folded.cos()
            elif folded > gamma/2:
                h = 1-r*(gamma-folded).cos()
            else:
                h = (r*folded.cos()).union(1-r*(gamma-folded).cos())
            low, high = enclose(h)
            peak_distance = min(enclose(folded)[0], enclose(2*gamma-folded)[0])
            trough_distance = enclose(abs(gamma-folded))[0]
            values.append((low, high, peak_distance, trough_distance))
        out = np.array(values, np.int64).T
        self.lattice[key] = out
        return out

    def translate(self, lo, hi):
        x0 = np.asarray(lo, dtype=np.int64)
        x1 = np.asarray(hi, dtype=np.int64)
        assert np.all((-UNIT < x0) & (x0 <= x1) & (x1 < UNIT))
        a, b = self.normal[:, :, 0], self.normal[:, :, 1]
        xx, yy = x0[:, None, :], x1[:, None, :]
        plus = np.maximum.reduce([xx*a, xx*b, yy*a, yy*b]).sum(axis=2)
        minus = np.minimum.reduce([xx*a, xx*b, yy*a, yy*b]).sum(axis=2)
        return minus//UNIT, -(-plus//UNIT)

    def boxes_supp(self, lo, hi, digital_scale):
        l = np.full((len(lo), self.m), UNIT//2, dtype=np.int64)
        h = l.copy()
        for n, i in BODY:
            left = [F(int(a), digital_scale) for a in lo[:, i]] if i >= 0 else [F(0)]*len(lo)
            right = [F(int(a), digital_scale) for a in hi[:, i]] if i >= 0 else left
            ll = np.stack([self.support(n, x)[0] for x in left])
            ul = np.stack([self.support(n, x)[1] for x in left])
            lr = np.stack([self.support(n, x)[0] for x in right])
            ur = np.stack([self.support(n, x)[1] for x in right])
            bmin, bmax = np.minimum(ll, lr), np.maximum(ul, ur)
            for j, (a, b) in enumerate(zip(left, right, strict=True)):
                if a == b:
                    continue
                _, _, peak, trough = self.support(n, (a+b)/2)
                w = (b-a)*UNIT/2
                bmin[j, trough <= int(w)] = enclose(1-self.radius[n])[0]
                bmax[j, peak <= int(w)] = enclose(self.radius[n])[1]
            assert digital_scale % UNIT == 0
            ratio = digital_scale//UNIT
            q0 = lo[:, i+1:i+3]//ratio
            q1 = -(-hi[:, i+1:i+3]//ratio)
            u0, u1 = self.translate(q0,q1)
            l = np.maximum(l, bmin+u0)
            h = np.maximum(h, bmax+u1)
        assert np.all(l <= h)
        return l.astype(np.int32), h.astype(np.int32)

    def specific_supp(self, p):
        low = np.full(self.m,UNIT//2,np.int64)
        high = low.copy()
        for n, i in BODY:
            rho = p[i] if i >= 0 else F(0)
            a, b, _, _ = self.support(n, rho)
            xy = p[i+1:i+3]
            xlo = np.array([[int((t*UNIT)//1) for t in xy]], np.int64)
            xhi = np.array([[ceil(t*UNIT) for t in xy]], np.int64)
            u, v = self.translate(xlo,xhi)
            low = np.maximum(low,a+u[0])
            high = np.maximum(high,b+v[0])
        return low, high

    def angular_bound(self, coordinate_max, scale):
        bounds = []
        for n, i in BODY:
            x, y = map(int, coordinate_max[i+1:i+3])
            translation_upper = F(sqrt_ceil(F((x*x+y*y)*UNIT**2,scale**2)),UNIT)
            derivative_upper = F(enclose(self.lipschitz[n])[1], UNIT)
            bounds.append(translation_upper + derivative_upper)
        return max(bounds)


def costs(outer_l, outer_u, center_l, center_u, perms):
    total, m = outer_l.shape
    assert center_u.shape == center_l.shape == (total, m)
    inf = np.iinfo(np.int32).max
    estimate = np.full((total,total),inf,np.int32)
    transport = np.zeros_like(estimate,dtype=np.uint8)
    for label, perm in enumerate(perms):
        for j, (l1,h1) in enumerate(zip(center_l[:,perm],center_u[:,perm],strict=True)):
            proposed = np.maximum((outer_u-l1).max(axis=1),(h1-outer_l).max(axis=1))
            improved = proposed < estimate[:,j]
            estimate[improved,j] = proposed[improved]
            transport[improved,j] = label
    return estimate,transport


def takes(estimate, gate, k):
    cells = estimate <= gate
    eligible = np.ones(len(estimate),bool)
    chosen = []
    for _ in range(k):
        hit = cells[eligible].sum(axis=0)
        j = int(hit.argmax())
        if not hit[j]:
            return None
        chosen.append(j)
        eligible &= ~cells[:,j]
        if not eligible.any():
            return chosen
    return None


def exact_two_cover(estimate, gate):
    bits = np.packbits((estimate <= gate).T, axis=1, bitorder='little')
    masks = [int.from_bytes(q.tobytes(), 'little') for q in bits]
    all_rows = (1 << len(estimate))-1
    for i, one in enumerate(masks):
        if one == all_rows:
            return [i]
        if not one:
            continue
        for j in range(i+1,len(masks)):
            if one|masks[j] == all_rows:
                return [i,j]
    return None


def shape_cover(estimate, k):
    # Only the attained cover is used. Greedy search need not give an optimum.
    low = int(estimate.min(axis=1).max())-1
    high = int(estimate.max(axis=0).min())
    search = lambda g: exact_two_cover(estimate,g) if k == 2 else takes(estimate,g,k)
    chosen = search(high)
    assert chosen is not None
    while high-low>1:
        gate = (high+low)//2
        trial = search(gate)
        if trial is None:
            low = gate
        else:
            chosen,high = trial,gate
    return chosen


def group_cover(estimate, groups):
    clients, operators = [], []
    for q in np.unique(groups):
        rows = np.flatnonzero(groups == q)
        j = rows[estimate[np.ix_(rows,rows)].max(axis=0).argmin()]
        clients.append(rows)
        operators.append(int(j))
    return clients,operators


def record(estimate,transport,tail, ids, selected, cells=None, shifts=None):
    if cells is None:
        assigned=np.argmin(estimate[:,selected],axis=1).astype(np.int32)
    else:
        assigned=np.zeros(len(estimate),np.int32)
        for j,group in enumerate(cells):
            assigned[group]=j
    reference = np.array(selected)[assigned]
    rows = np.arange(len(estimate))
    radius=estimate[rows,reference].astype(np.int64)+tail
    symmetry=transport[rows,reference]
    cluster_radii = [int(radius[assigned == j].max(initial=0)) for j in range(len(selected))]
    data = dict(centers=[ids[j] for j in selected], center_radius=cluster_radii,
                global_radius=int(radius.max()),
                rows=[[int(q),int(e),int(v)] for q,e,v in zip(assigned,symmetry,radius,strict=True)])
    if shifts is not None:
        data['translation_numerators'] = shifts[rows,reference].tolist()
    return data


def translated_costs(s, l, u, c, d, matrix, move, tail, candidates):
    # Numerical first moments only propose a rigid shift. Its entire certificate
    # uses interval normal components and integer projections, including the tail.
    point = (s.normal[:,:,0]+s.normal[:,:,1])//2
    assert np.all(s.normal[:,:,1]-s.normal[:,:,0] == 2)
    template_center, leaf_center = (c+d)/2, (l+u)/2
    out = matrix+tail
    shifts = np.zeros((*out.shape,2),np.int32)
    for j in candidates:
        perm=s.perm[move[:,j]]
        med=template_center[j,perm]
        center_proposal=2*(leaf_center-med)@(point/UNIT)/s.m
        base_l, base_u=c[j,perm],d[j,perm]
        for step in (1,2,3,4,5,6):
            t=np.rint(center_proposal*step/4).astype(np.int64)
            size=np.abs(t).sum(axis=1)
            assert np.all(size < UNIT//4)
            projection=t[:,0:1]*point[:,0]+t[:,1:2]*point[:,1]
            ts0=(projection-size[:,None])//UNIT
            ts1=-(-(projection+size[:,None])//UNIT)
            distance=np.maximum((u-base_l-ts0).max(axis=1),(base_u+ts1-l).max(axis=1))
            correction=-(-22*size//(7*s.m))
            value=distance+tail+correction
            better=value<out[:,j]
            out[better,j]=value[better]
            shifts[better,j]=t[better]
    return out,shifts


def escape_check(s, main):
    from certify_five_harmonic import placement_area
    b, bhash = checked(main/'five-body-contact-basin.json',
                      '282fafa2a4d37f72871c082b9e3e3e71c19fdb9ee4889cfa1a08fef19dcd150a')
    p, phash = checked(main/'competing-fits-input.json',
                       '31aaf17ffb9e33b4a50050614c37c78fdcfd9976360d561bd42eb1dd26070f96')
    r, rhash = checked(main/'competing-fits-replay.json',
                       'b2ec3077aa04306c05a0f25829ddf6b0ca6f6e049f60776ed659e9158dfcacca')
    binput, bihash = checked(main/'basin-083945-input.json',
                              'a1d5f8b0889ea0fca0dcf5a1566648ed5e27c63f5783df9e385c3ba4782c38f6')
    bcert, bcert_hash = checked(main/'basin-083945-parent.json',
                               'bf51d740089e047540ac266541a505594f71fef8c894eb6a970adf79f2569537')
    assert binput['center'] == b['raw_center'] and binput['terms'] == p['terms']
    assert F(bcert['lower']) > F('.83945005239')
    f = r['five_with_s_strictly_inside']
    upper = F(f['area_upper'])
    assert upper < F('.837657369') and F(f['minimum_support_slack']) > 0
    low_sample = F('.000534')
    assert (1+4*low_sample)*upper + PI_UP*low_sample**2 < F('.83945005239')
    center = list(map(F,b['raw_center']))
    esc = list(map(F,p['zero_gain_five']['params']))
    esc_lower = F(placement_area(esc, s.pi, k=128)['lower'])
    contact_upper = F(placement_area(center, s.pi, k=128)['upper'])
    assert esc_lower > F('.837295006222')
    assert contact_upper < F('.837657369')
    l1,h1 = s.specific_supp(center)
    l2,h2 = s.specific_supp(esc)
    lo = max(int((l2-h1).max()), int((l1-h2).max()))
    hi = max(int((h2-l1).max()), int((h1-l2).max()))
    mat = np.array([[ceil(max(abs(c),abs(e))*2**60) for c,e in zip(center,esc,strict=True)]])
    lip = s.angular_bound(mat[0], 2**60)
    tail = ceil(2*lip*PI_UP*UNIT/s.m)
    direction = (s.index+s.m//2)%s.m
    wl, wh, el, eh = l1+l1[direction], h1+h1[direction], l2+l2[direction],h2+h2[direction]
    difference = 2*UNIT
    for sign in (1,-1):
        for j in range(s.m):
            rot=(sign*s.index+j)%s.m
            sample=max(int((wl-eh[rot]).max()),int((el[rot]-wh).max()))
            difference=min(difference,sample)
    width_distance_floor = (difference-tail)//2
    assert width_distance_floor > low_sample*UNIT
    return dict(same_frame_lower=lo,same_frame_upper=hi+tail, arbitrary_isometry_lower=width_distance_floor,
                area_only_arbitrary_isometry_lower=str(low_sample),area_lower=str(esc_lower),
                area_upper=str(upper), contact_area_upper=str(contact_upper),
                input_sha256={'contact':bhash,'escape':phash,'escape_receipt':rhash,
                              'contact_six_input':bihash,'contact_six_replay':bcert_hash})


def transfer_budget(s, target):
    r, radius = F(1,2), F(703,1000)
    # Tangent triangle outside the radius-1/2 disk minus the overlaid sector.
    d = exact(radius*radius-r*r).sqrt()
    area = s.pi/4+exact(r)*d-(exact(r/radius)).acos()/4
    assert area > exact(F('.83732'))
    factor = F(sqrt_ceil((4*radius**2-1)*UNIT**2), UNIT)
    return sqrt_ceil((target-GLOBAL_LOW)*factor*UNIT**2), radius


def full_box_obstruction(s, lo, hi, ids):
    # Explicit corner placements in one raw leaf. No area condition is imposed.
    a,b=2*lo,2*hi
    a[:,[2,5,8]]=b[:,[2,5,8]]=lo[:,[2,5,8]]+hi[:,[2,5,8]]
    h0,h1=s.boxes_supp(a,b,2*RAW)
    i,j=np.unravel_index(np.argmax(h1-h0),h0.shape)
    positive=s.normal[j,:,0] >= 0
    upper=a[i].copy()
    lower=b[i].copy()
    for _, k in BODY:
        upper[k+1:k+3]=np.where(positive,b[i,k+1:k+3],a[i,k+1:k+3])
        lower[k+1:k+3]=np.where(positive,a[i,k+1:k+3],b[i,k+1:k+3])
    assert np.all(a[i]<=upper) and np.all(upper<=b[i])
    assert np.all(a[i]<=lower) and np.all(lower<=b[i])
    l,_=s.specific_supp([F(int(x),2*RAW) for x in upper])
    _,h=s.specific_supp([F(int(x),2*RAW) for x in lower])
    separation=int(l[j]-h[j])
    assert separation > F('.02')*UNIT
    return dict(raw_leaf=ids[i],normal_index=int(j),corner_poses=[upper.tolist(),lower.tolist()],
                separation_lower=separation,single_center_radius_lower=separation//2)


def verify(a):
    ctx.prec=a.bits
    ids,lo,hi,group,provenance=read_leaves(a.computation,a.tree)
    assert len(ids)==501 and np.all(lo<hi)
    s=Samples(a.normals)
    x1,x2=s.boxes_supp(lo,hi,RAW)
    c1,c2=s.boxes_supp(lo+hi,lo+hi,2*RAW)
    assert np.all(x1<=c2) and np.all(c1<=x2)
    norm_max=np.maximum(np.abs(lo),np.abs(hi)).max(axis=0)
    lip=s.angular_bound(norm_max, RAW)
    tail=ceil(2*lip*PI_UP*UNIT/s.m)
    m1,t1=costs(x1,x2,c1,c2,s.perm[:1])
    m2,t2=costs(x1,x2,c1,c2,s.perm)
    assert np.all(m2<=m1)
    methods={}
    def use(name,matrix,moves,k):
        methods[name]=record(matrix,moves,tail,ids,shape_cover(matrix,k))
    for k in (1,2,3,4,6,8,12,16,24,32):
        use(f'nearest_D3_{k}',m2,t2,k)
        use(f'nearest_fixed_{k}',m1,t1,k)
    cells,selector=group_cover(m2,group)
    methods['R138_components_3']=record(m2,t2,tail,ids,selector,cells)
    cut=3684498829194
    for m in (1,2,3,4,6,8):
        index=np.minimum(m-1,(lo[:,2]+hi[:,2])*m//(2*cut))
        label=(group!=0)*m+index
        cells,selector=group_cover(m2,label)
        methods[f'angle_histogram_{m}']=record(m2,t2,tail,ids,selector,cells)
    candidates={tuple(key) for kind,r in methods.items()
                if kind in ('nearest_D3_1','nearest_D3_2','nearest_D3_3','nearest_D3_4',
                            'nearest_D3_6','nearest_D3_8','nearest_D3_12','R138_components_3')
                for key in r['centers']}
    selected=[j for j,key in enumerate(ids) if tuple(key) in candidates]
    m3,shifts=translated_costs(s,x1,x2,c1,c2,m2,t2,tail,selected)
    assert np.all(m3<=m2+tail)
    for k in (1,2,3,4,6,8,12):
        methods[f'nearest_translated_{k}']=record(m3,t2,0,ids,shape_cover(m3,k),shifts=shifts)
    supports={tuple(key) for r in methods.values() for key in r['centers']}
    candidates=[{ 'leaf':ids[j], 'pose_numerator':[int(q) for q in (lo[j]+hi[j])] }
                for j in range(len(lo)) if tuple(ids[j]) in supports]
    target=F('.837295006222')
    transfer,radius=transfer_budget(s,target)
    if a.emit:
        np.savez_compressed(a.emit,ids=np.array(ids),lo=lo,hi=hi,group=group,box_lo=x1,box_hi=x2,
                           midpoint_lo=c1,midpoint_hi=c2,fixed=m1,dihedral=m2,dihedral_move=t2)
    return dict(area_target=str(target),lower_area=str(GLOBAL_LOW),coordinate_den=RAW,prototype_den=2*RAW,
                support_den=UNIT,normal_count=s.m,angular_tail=tail,angular_lipschitz=str(lip),
                larger_host_extra=transfer,host_circumradius=str(radius),methods=methods,prototypes=candidates,leaves=ids,
                raw_box_sha256=digest(lo.tobytes()+hi.tobytes()),digest={**provenance,'source':digest(Path(__file__).read_bytes())},
                escape=escape_check(s,a.computation),
                fixed_frame_whole_leaf_obstruction=full_box_obstruction(s,lo,hi,ids))


if __name__=='__main__':
    p=ArgumentParser(description=__doc__)
    p.add_argument('computation',type=Path)
    p.add_argument('tree',type=Path)
    p.add_argument('--normals',type=int,default=3072)
    p.add_argument('--bits',type=int,default=160)
    p.add_argument('--emit',type=Path)
    print(json.dumps(verify(p.parse_args()),separators=(',',':')))
