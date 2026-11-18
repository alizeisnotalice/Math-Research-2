"""Exact checks of the parity and disjoint-block inequalities for n <= 4."""
import json
from math import factorial
import platform


def bits(mask):
    while mask:
        bit = mask & -mask
        yield bit
        mask ^= bit


def index(bit):
    return bit.bit_length() - 1


def case(n):
    N = n*n
    adj = [sum(1 << v for v in range(N)
               if abs(u%n-v%n)+abs(u//n-v//n)==1) for u in range(N)]
    masks = 1 << N
    all_vertices = masks - 1
    even = sum(1 << v for v in range(N) if (v%n+v//n)%2 == 0)
    odd = all_vertices ^ even
    around = [0]*masks
    moments = [0]*masks
    signed = [0]*masks
    upper_num = [1]*masks
    upper_den = [1]*masks
    moments[0] = signed[0] = factorial(N)
    for mask in range(1,masks):
        v_bit = mask & -mask
        v = index(v_bit)
        rest = mask ^ v_bit
        around[mask] = around[rest] | adj[v]
        degree = adj[v].bit_count()
        upper_num[mask] = upper_num[rest]*degree
        upper_den[mask] = upper_den[rest]*(degree+4)
        if mask & around[mask]:
            continue
        numerator = sum(moments[mask^b] for b in bits(mask))
        value, remainder = divmod(numerator,(mask|around[mask]).bit_count())
        assert remainder == 0
        moments[mask] = value
        signed[mask] = -value if mask.bit_count()%2 else value
    C = signed.copy()
    for v in range(N):
        v_bit = 1 << v
        for mask in range(v_bit,masks):
            if mask & v_bit:
                C[mask] += C[mask ^ v_bit]
    fact = factorial(N)
    suffix = [0]*masks
    suffix[-1] = 1
    for S in range(all_vertices-1,0,-1):
        suffix[S] = sum(suffix[S|v] for v in bits(around[S]&~S))

    single_extensions = components_checked = 0
    parity_strict = 0
    for U in range(masks):
        assert 0 <= C[U] <= fact
        e, o = U&even, U&odd
        assert C[U]*fact <= C[e]*C[o]
        parity_strict += C[U]*fact < C[e]*C[o]
        if not e & around[o]:
            assert C[U]*fact == C[e]*C[o]
        assert C[U]**4*upper_den[U] <= fact**4*upper_num[U]
        for bit in bits(all_vertices & ~(U|around[U])):
            d = adj[index(bit)].bit_count()
            assert (d+1)*C[U|bit] >= d*C[U]
            single_extensions += 1
        remaining = U
        product = 1
        count = 0
        while remaining:
            component = remaining & -remaining
            next_set = component | (around[component]&remaining)
            while next_set != component:
                component = next_set
                next_set = component | (around[component]&remaining)
            count += 1
            product *= C[component]
            remaining ^= component
        if count > 1:
            assert C[U]*fact**(count-1) >= product
            components_checked += 1
        if U:
            # U also serves as the independent crosscheck's exposed set.
            forbidden = all_vertices & ~(U|around[U])
            assert suffix[U]*fact == C[forbidden]*factorial(N-U.bit_count())
    return dict(n=n, tested_masks=masks, parity_strict=parity_strict,
                single_extensions=single_extensions, components_checked=components_checked,
                suffix_masks=masks-1)


if __name__ == '__main__':
    print(json.dumps(dict(python=platform.python_version(), checks=[case(n) for n in (2,3,4)]),indent=2))
