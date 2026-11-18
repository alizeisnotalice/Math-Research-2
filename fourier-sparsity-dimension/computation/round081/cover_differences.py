"""Check the weighted uncovered-mass inequality over F2^3 by integer Walsh sums."""
import itertools
import json
import platform

def spaces(n):
    result={frozenset({0})}
    for v in range(1,1<<n):
        result.update([s|{x^v for x in s} for s in result])
    return sorted((sorted(s) for s in result),key=lambda s:(len(s),s))

if __name__=='__main__':
    n=3
    size=1<<n
    masks=spaces(n)
    spectra=set()
    for truth in range(1<<size):
        values=[1-2*((truth>>x)&1) for x in range(size)]
        spectra.add(tuple(abs(sum(v*(-1)**((x&e).bit_count()%2)
                                  for x,v in enumerate(values))) for e in range(size)))
    counts=[]
    for m in range(1,5):
        checked=tight=0
        for cover in itertools.combinations_with_replacement(masks,m):
            outside=set(range(size)).difference(*cover)
            for e in outside:
                weights=[2**sum((e^g) in u for u in cover) if g in outside else 0
                         for g in range(size)]
                for spec in spectra:
                    if spec[e]:
                        total=sum(w*a for w,a in zip(weights,spec))
                        assert total>=2*size
                        checked+=1
                        tight+=total==2*size
        counts.append(dict(m=m,checked=checked,tight=tight))
    print(json.dumps(dict(python=platform.python_version(),n=n,
                     spectra=len(spectra),spaces=len(masks),counts=counts),indent=2))
