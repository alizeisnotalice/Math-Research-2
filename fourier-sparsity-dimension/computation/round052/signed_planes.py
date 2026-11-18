"""Exact signed-band and crossover checks for tilted two-plane partitions."""
from fractions import Fraction
import json


def parity(x):
    return x.bit_count() & 1


def fwht(values):
    out = list(values)
    step = 1
    while step < len(out):
        for a in range(0, len(out), 2*step):
            for j in range(a, a+step):
                v, w = out[j], out[j+step]
                out[j], out[j+step] = v+w, v-w
        step *= 2
    return out


def quad(x, bits):
    return sum(((x >> i) & 3) == 3 for i in range(0, bits, 2)) % 2


def check(k):
    normal = (1, 2, 4, 7)
    selector_count = 32 * 4**k
    sign = []
    cells = {}
    for addr in range(selector_count):
        x, w, u = addr & 7, (addr >> 3) & 3, addr >> 5
        xi = [(x >> i) & 1 for i in range(3)]
        q = xi[0]*xi[1] + xi[0]*xi[2] + xi[1]*xi[2] + quad(w, 2) + quad(u, 2*k)
        sign.append((-1)**q)
        cells.setdefault((u, w, parity(normal[w] & x)), []).append(addr)
    m = len(cells)
    weighted_signs = [0]*selector_count
    for points in cells.values():
        origin, b, c, _ = points
        d1, d2 = origin ^ b, origin ^ c
        vals = [sign[origin], sign[origin ^ d1], sign[origin ^ d2], sign[origin ^ d1 ^ d2]]
        local = fwht(vals)
        assert all(abs(x) == 2 for x in local)
        for t in range(selector_count):
            eta = parity(t & d1) + 2*parity(t & d2)
            weighted_signs[t] += (-1 if parity(t & origin) else 1) * (local[eta] // 2)
    global_spectrum = fwht(sign)
    assert all(a*2*m == b*selector_count for a, b in zip(global_spectrum, weighted_signs))
    s_band = sum(x != 0 for x in global_spectrum)
    signed_bound = Fraction((2*m)**2, max(map(abs,weighted_signs))**2)
    assert s_band == signed_bound == 16*4**k
    # Tensor with the two-data-bit bent sign (-1)^(y0*y1).
    joint = fwht([s*(-1 if y == 3 else 1) for y in range(4) for s in sign])
    assert sum(x != 0 for x in joint) == 4*s_band
    chosen_rank_bound = 5+2*k+2+2
    local_D = 6
    refund = m*local_D - chosen_rank_bound
    unsigned_needed = (m-1)*local_D
    signed_needed = m*local_D - 3*(4*2**k)
    assert refund < unsigned_needed and refund >= signed_needed
    return dict(k=k, cells=m, shared_bands=4, lifts_per_active_cell=4,
                unsigned_lower_bound=4, signed_lower_bound=int(signed_bound),
                actual_global_sparsity=4*s_band, rank_refund=refund,
                unsigned_refund_needed=unsigned_needed, signed_refund_needed=signed_needed)


if __name__ == '__main__':
    print(json.dumps([check(k) for k in range(4)], indent=2))
