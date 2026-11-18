"""Exhaust the four-bit truth tables with integer autocorrelation and Walsh spectra."""
from collections import Counter
import json
import platform

def walsh(row):
    row = row.copy()
    n = 1
    while n < len(row):
        for i in range(0,len(row),2*n):
            for j in range(i,i+n):
                a,b = row[j],row[j+n]
                row[j],row[j+n] = a+b,a-b
        n *= 2
    return row

def expand(mask, label):
    return mask | sum(1 << (i^label) for i in range(16) if mask >> i & 1)

def main():
    hist = Counter()
    for truth in range(1 << 16):
        row = [1-2*((truth>>i)&1) for i in range(16)]
        spectrum = walsh(row)
        labels = [i for i,v in enumerate(spectrum) if v and i]
        corr = [sum(row[i]*row[i^j] for i in range(16)) for j in range(16)]
        assert walsh(corr) == [x*x for x in spectrum]
        k = len(set(corr))-1
        levels = [1] * (k+1)
        span = 1
        for e in labels:
            span = expand(span,e)
            for j in range(k,0,-1):
                levels[j] |= expand(levels[j-1],e)
        assert levels[k] == span
        hist[k+1] += 1
    return dict(python=platform.python_version(),tables=1<<16,
                by_autocorrelation_value_count=dict(sorted(hist.items())))

if __name__ == '__main__':
    print(json.dumps(main(),indent=2))
