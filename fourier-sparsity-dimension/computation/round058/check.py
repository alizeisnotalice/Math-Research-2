import json
from itertools import combinations
def basis(xs):
    rows={}
    for x in xs:
        while x and x.bit_length() in rows: x ^= rows[x.bit_length()]
        if x:rows[x.bit_length()]=x
    return rows
def rank(xs):return len(basis(xs))
def W(a):
    a=a.copy(); step=1
    while step<len(a):
        for p in range(0,len(a),2*step):
            for i in range(p,p+step):
                a[i],a[i+step]=a[i]+a[i+step],a[i]-a[i+step]
        step*=2
    return a
q0=[int((x&3)==3)^int((x&12)==12) for x in range(16)]
lift=[int((x&5)==5) for x in range(16)]
for y in (0,1):
    assert set(map(abs,W([1-2*(a^(y*b)) for a,b in zip(q0,lift)])))=={4}
rels=[sum(1<<v for v in quad) for quad in combinations(range(16),4)
      if quad[0]^quad[1]^quad[2]^quad[3]==0]
r=rank(rels)
A=sum(val<<i for i,val in enumerate(lift))
nonzero=any((rel&A).bit_count()%2 for rel in rels)
true=rank([x^(1<<4)^(lift[x]<<5) for x in range(16)])
assert (r,nonzero,true)==(11,True,6)
# Two packets of two basis labels, each with its own selector, and one outer bit.
signals=[]
for x in range(128):
    g=x&1
    u=(x>>(1+3*g))&1
    bit=(x>>(2+3*g+u))&1
    signals.append(1-2*bit)
glob=[i for i,c in enumerate(W(signals)) if c]
B={x>>1 for x in glob}
assert len(glob)==16 and rank(B)==6
print(json.dumps(dict(expected_outer_intersection_dimension=r-int(nonzero),chosen_relations_rank=r,
        full_walsh_rank=true,two_packet_examples=dict(s=len(glob),data_rank=rank(B),
        full_rank=rank(glob))),indent=2))
