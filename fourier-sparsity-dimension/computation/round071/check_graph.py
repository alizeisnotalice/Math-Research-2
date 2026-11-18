"""Direct integer truth tables, Walsh coefficients and graph differences."""
from collections import Counter
import json
import platform


def run(m, a, polynomial):
    n = 1 << m
    def mul(x, y):
        value = 0
        while y:
            if y & 1:
                value ^= x
            y >>= 1
            x <<= 1
            if x & n:
                x ^= polynomial
        return value
    def trace(x):
        value = 0
        for _ in range(m):
            value ^= x
            x = mul(x,x)
        assert value in (0,1)
        return value
    def dual(x):
        return sum(trace(mul(x,1<<i)) << i for i in range(m))

    row = []
    for t in range(1<<a):
        c = (1<<a) ^ t
        inv = next(y for y in range(n) if mul(c,y)==1)
        for x in range(n):
            for y in range(n):
                row.append(1 - 2*trace(mul(mul(x,y),inv)))
    width = 1
    while width < len(row):
        for i in range(0,len(row),2*width):
            for j in range(i,i+width):
                l,r = row[j],row[j+width]
                row[j],row[j+width] = l+r,l-r
        width *= 2
    expected = {}
    for p in range(n):
        for q in range(n):
            product = mul(p,q)
            head = dual(product) & ((1<<a)-1)
            index = (head << (2*m)) | (dual(p) << m) | dual(q)
            expected[index] = (1<<(m+a)) * (1-2*trace(mul(1<<a,product)))
    assert {i:v for i,v in enumerate(row) if v} == expected
    counts = Counter(x^y for x in expected for y in expected)
    for w in range(1,len(row)):
        assert counts[w] == (n*n//(1<<a) if w & (n*n-1) else 0)
    return dict(m=m, a=a, field_polynomial=polynomial, truth_entries=len(row),
                spectrum=len(expected), good_differences=len(counts)-1,
                nonzero_intersection_size=n*n//(1<<a))


if __name__ == '__main__':
    print(json.dumps(dict(python=platform.python_version(),
                          cases=[run(3,2,0b1011),run(4,2,0b10011),run(4,3,0b10011)]),
                     indent=2))
