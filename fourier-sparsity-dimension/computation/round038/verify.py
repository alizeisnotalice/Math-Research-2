"""Check the 108-way partition, counts and one minimum witness per shard."""
from pathlib import Path
from itertools import product

here = Path(__file__).resolve().parent
rows = [line.split() for name in ("mini", "macbook")
        for line in (here / f"{name}-complete.txt").read_text().splitlines()]
assert sorted(int(r[0]) for r in rows) == list(range(108))
assert all(int(r[1]) == 108 for r in rows)
def mask(fn):
    return sum(int(fn([(x >> i) & 1 for i in range(6)])) << x for x in range(64))
base = mask(lambda x: x[0]*x[1] ^ x[2]*x[3] ^ x[4]*x[5])
def switch(x, t):
    a,b,c,d,e = t
    return (1^x[5])*((1^x[3])*(a^b*x[0]^c*x[1]) ^ x[3]*(d*x[0]^e*x[1]))
words = [mask(lambda x: switch(x,t)) for t in product(range(2),repeat=5)]
assert len(set(words))==32
support = 0
for word in words:
    support |= word
assert support.bit_count()==28
chars = [sum(((i & x).bit_count()%2)<<x for x in range(64)) for i in range(64)]
for _,_,_,count,size,hexword in rows:
    if int(count):
        v = int(hexword,16)
        assert v not in words and (v|support).bit_count()==int(size)
        assert all((v^word^base^c).bit_count() in (28,36) for word in words for c in chars)
total=sum(int(r[2]) for r in rows)
independent=sum(int(r[3]) for r in rows)
best=min(int(r[4]) for r in rows)
assert (total,independent,best)==(688128,688096,30)
print(f"shards={len(rows)}, accepted={total}, independent={independent}, minimum={best}")
