"""R003: independent full-law calculation and rational exit-extension checks."""
from fractions import Fraction as F
from math import prod
from platform import python_version

def convolution(row):
    sums = {}
    for i, x in row.items():
        for j, y in row.items():
            sums[i+j] = sums.get(i+j, 0) + x*y
    return sums

def run(law, depth):
    laws = [law]
    for _ in range(depth):
        sums = convolution(laws[-1])
        new = {k-1:v for k,v in sums.items() if k>1}
        new[0] = sums.get(0,0)+sums.get(1,0)
        laws.append(new)
    a = [{k:2**k*p for k,p in row.items() if k>0} for row in laws]
    sources = [convolution(row) for row in a]
    def weight(u,v):
        return prod(row[0] for row in laws[u:v])
    checks = lost = shift_checks = 0
    for u in range(depth):
        for v in range(u+1,depth+1):
            # Mass that arrived at zero at i+1 is kept at i+1-v at time v.
            ext = dict(a[v])
            for i in range(u,v):
                ext[i+1-v] = weight(i,v)*a[i].get(1,0)
            kernel = convolution(ext)
            for j in range(u,v):
                for x,ax in a[j].items():
                    assert ext.get(x+j-v,0) >= weight(j,v)*ax
                    shift_checks += 1
            for n in range(v,depth+1):
                for k in range(1,7):
                    actual = sum(weight(j+1,n)*sources[j].get(n-j+k,0) for j in range(u,v))
                    slab = range(n-2*v+u+k,n-v+k)
                    bound = sum(kernel.get(s,0) for s in slab)/weight(u,v)**2
                    assert actual <= bound
                    checks += 1
                    poskernel = convolution(a[v]) if k==1 and n==v else None
                    if poskernel is not None:
                        lost += int(sum(poskernel.get(s,0) for s in slab)/weight(u,v)**2 < actual)
    return checks,shift_checks,lost

if __name__ == '__main__':
    print(python_version())
    for law in ({0:F(4,5),2:F(1,5)}, {0:F(8,17),1:F(1,2),3:F(1,34)}):
        print(run(law,6))
