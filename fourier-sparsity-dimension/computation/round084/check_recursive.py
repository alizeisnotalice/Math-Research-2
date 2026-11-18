"""python3 -B check_recursive.py. Independent integer checks of the selector budgets."""
import json
import platform
from fractions import Fraction


def rank(columns):
    basis = {}
    for v in columns:
        while v:
            key = v.bit_length()
            if key not in basis:
                basis[key] = v
                break
            v ^= basis[key]
    return len(basis)


def D(support):
    return 1 + rank(1 + 2*v for v in support)


def span(axes):
    labels = [0]
    for bit in axes:
        labels += [x ^ bit for x in labels]
    return labels


class Space:
    def __init__(self):
        self.count = 0

    def axes(self, n):
        bits = [1 << (self.count + j) for j in range(n)]
        self.count += n
        return bits


def paths(space, n):
    nodes = space.axes(n - 1)
    leaves = []
    for index in range(n, 2*n):
        path = []
        while index > 1:
            index //= 2
            path.append(nodes[index-1])
        leaves.append(path)
    return leaves


def masks(space, M):
    j = M.bit_length() - 1
    T = M*M
    t, k = (2*T + 2)//5, (M+3)//7
    outer, coarse, fine = space.axes(2), space.axes(2), space.axes(6*j)
    pad = set()
    for i, path in enumerate(paths(space, M)):
        guard, = space.axes(1)
        pad.update(guard ^ x for x in span(path + space.axes(3 if i < k else 2)))
    K = 4*M*(M+k)
    assert len(pad) == K
    output = []
    prescribed = []
    for path in paths(space, t):
        guard, = space.axes(1)
        output.append({guard ^ x ^ y for x in pad for y in span(outer + path)})
        prescribed.append(4*K*2**len(path))
    for incident in ((1, 2), (1, 3), (2, 3)):
        # Complement of two distinct coarse points is a line.
        normal = sum(bit for r, bit in enumerate(coarse)
                     if (incident[0] ^ incident[1]) & (1 << (1-r)))
        for path in paths(space, T):
            guard, = space.axes(1)
            output.append({guard ^ y for y in span(outer + [normal] + path)})
            prescribed.append(8*T)
    combined = set()
    for edge in ((1, 2), (1, 3), (2, 3)):
        guard, = space.axes(1)
        normal = sum(bit for r, bit in enumerate(outer)
                     if (edge[0] ^ edge[1]) & (1 << (1-r)))
        combined.update(guard ^ y for y in span([normal] + coarse + fine))
    output.append(combined)
    prescribed.append(24*M**6)
    assert list(map(len, output)) == prescribed
    return output


def seed(space, family):
    if family == 'constant':
        return {0}
    if family == 'character':
        return set(space.axes(1))
    if family == 'bent':
        return set(span(space.axes(2)))
    if family == 'full':
        return set(span(space.axes(3)))
    a, b, x = space.axes(3)
    return {a, b, a ^ x, b ^ x}


def params(M):
    j = M.bit_length()-1
    T = M*M
    t, k = (2*T+2)//5, (M+3)//7
    K = 4*M*(M+k)
    e = 4*M+k+6*j
    rate = Fraction(3,8) + Fraction(1,4*K) + Fraction((2+e)**2,96*M**6)
    return T, t, K, e, rate


def test_spectrum(M, pattern, bg):
    space = Space()
    T, t, K, e, rate = params(M)
    bands = masks(space, M)
    kinds = [pattern[i % len(pattern)] for i in range(len(bands)-1)] + [bg]
    children = [seed(space, kind) for kind in kinds]
    full = set()
    total = sum(len(xs)*len(ys) for xs, ys in zip(bands, children))
    for band, child in zip(bands, children):
        full.update(x ^ y for x in band for y in child)
    assert len(full) == total
    Ds = list(map(D,children))
    parent = D(full)
    assert parent == sum(Ds) + e
    relative = Fraction(parent*parent, total) / max(
        Fraction(dd*dd,len(c)) for dd,c in zip(Ds,children))
    assert relative <= rate <= Fraction(83,192)
    return dict(M=M, independent_slots=len(bands), pattern=pattern, background=bg,
                actual_s=total, actual_D=parent, squared_relative=str(relative))


def limit_test():
    rows = []
    for j in [1,2,3,4,6,9,16,32]:
        M = 2**j
        T,t,K,e,upper = params(M)
        assert Fraction(3,8) < upper <= Fraction(83,192)
        # Full heap leaves have depths log2(lower) and log2(lower)+1.
        lower = 1 << (t.bit_length()-1)
        sA_const = 4*K*lower*(3*t-2*lower)
        L = 2**(5*j+28)
        parent = 6*T*L + 2*(t+1) + e
        root_s = 24*T*T*L*L+sA_const+24*M**6
        ratio = Fraction(parent*parent,4*root_s)
        assert abs(ratio-Fraction(3,8)) < Fraction(1,10**6)
        assert ratio <= upper
        rows.append(dict(M=M, rho_squared=str(upper), rational_witness=str(ratio)))
    return rows


if __name__ == '__main__':
    print(json.dumps(dict(python=platform.python_version(), seed=None,
                          tolerance=0, finite_fixtures=[
                              test_spectrum(2,['constant','character','bent','full','tree'],'tree'),
                              test_spectrum(2,['tree','full'],'bent'),
                              test_spectrum(4,['tree','bent','constant'],'character')],
                          parameter_checks=limit_test()), indent=2))
