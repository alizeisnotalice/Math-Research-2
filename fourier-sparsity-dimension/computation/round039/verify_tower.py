"""Independently verify seed ranks and Walsh magnitudes for the nested witnesses."""
from pathlib import Path

def truth(fn):
    return sum(fn([(x>>i)&1 for i in range(6)]) << x for x in range(64))
base=truth(lambda x: x[0]*x[1]^x[2]*x[3]^x[4]*x[5])
table=[truth(f) for f in (
    lambda x:(1^x[5])*(1^x[3]),
    lambda x:(1^x[5])*(1^x[3])*x[0],
    lambda x:(1^x[5])*(1^x[3])*x[1],
    lambda x:(1^x[5])*x[3]*x[0],
    lambda x:(1^x[5])*x[3]*x[1])]
table.append(0x3c5a3c5a)
chars=[sum(((x&i).bit_count()%2)<<x for x in range(64)) for i in range(64)]
phases=[base]
support=0
for w in table:
    support |= w
    phases += [p^w for p in phases]
here=Path(__file__).resolve().parent
for d in range(7,11):
    rows=[line.split() for n in (12,15) for line in (here / f'tower-{d}-{n}.txt').read_text().splitlines()]
    assert sorted(int(r[0]) for r in rows)==list(range(108))
    eligible=[(int(r[4]),int(r[5],16)) for r in rows if int(r[3])]
    for size,w in eligible:
        assert (support|w).bit_count()==size and base^w not in phases
        assert all((p^w^c).bit_count() in (28,36) for p in phases for c in chars)
    size,w=min(eligible)
    support |= w
    phases += [p^w for p in phases]
    assert len(set(phases))==2**d
    print(d,size,format(w,'016x'),sum(int(r[3]) for r in rows))
