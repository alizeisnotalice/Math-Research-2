"""Round 2. Exact finite-orbit replay of the transport and frozen-block inequalities."""
from fractions import Fraction as Q
from math import prod
from platform import python_version

def next_law(p):
    out = [Q(0)] * max(1, 2*len(p)-2)
    for i, x in enumerate(p):
        for j, y in enumerate(p):
            out[max(0, i+j-1)] += x*y
    return out

def value(a, i):
    return a[i] if 0 < i < len(a) else Q(0)

def audit(law, last):
    p = [law]
    for _ in range(last):
        p.append(next_law(p[-1]))
    a = [[Q(0)] + [2**k*v for k, v in enumerate(row[1:], 1)] for row in p]
    source = []
    for row in a[:-1]:
        c = [Q(0)] * (2*len(row)-1)
        for i in range(1, len(row)):
            for j in range(1, len(row)):
                c[i+j] += row[i]*row[j]/2
        source.append(c)
    checks = nonzero = identities = 0
    for n in range(1,last+1):
        for k in range(1,7):
            w = lambda u,v: prod(p[i][0] for i in range(u,v))
            terms = [w(j+1,n)*value(source[j], n-j+k) for j in range(n)]
            assert value(a[n],k) == w(0,n)*value(a[0],n+k)+sum(terms)
            identities += 1
            for u in range(n):
                for v in range(u+1,n+1):
                    h = v-u
                    low = n-u+k
                    # Count each pair once in the frozen slab.
                    frozen = sum((a[u][x]*a[u][y]/2
                                  for x in range(h+1,len(a[u]))
                                  for y in range(h+1,len(a[u]))
                                  if low <= x+y < low+h), Q(0))
                    bound = w(u,n)**2*frozen
                    assert sum(terms[u:v]) >= bound
                    checks += 1
                    nonzero += int(bound>0)
    return identities,checks,nonzero

if __name__ == '__main__':
    print(python_version())
    for law in ([Q(4,5),Q(0),Q(1,5)], [Q(8,17),Q(1,2),Q(0),Q(1,34)]):
        assert sum(law)==1 and sum((k-1)*2**k*x for k,x in enumerate(law))==0
        print(audit(law,6))
