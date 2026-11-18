"""Private-pair grid selectors, raw integer FWHT, and a weighted four-leaf certificate."""
from fractions import Fraction
from itertools import permutations
import json
import platform

def walsh(row):
    a = list(row)
    stride = 1
    while stride < len(a):
        for start in range(0,len(a),2*stride):
            for j in range(start,start+stride):
                x,y=a[j],a[j+stride]
                a[j],a[j+stride]=x+y,x-y
        stride*=2
    return a

def grid(m):
    t=m*m
    # Each tag is a vector in the independent e,f,u coordinate system.
    bands={}
    for i in range(m):
        for j in range(m):
            index=i*m+j
            e,f,u=1<<i,1<<(m+j),1<<(2*m+index)
            for beta,sign in ((e,1),(f,1),(u,1),(e^f^u,-1)):
                bands.setdefault(beta,[0]*t)[index]=sign
    s=sum(sum(v!=0 for v in walsh(row)) for row in bands.values())
    assert s==2*t*t+2*m*m
    if m==2:
        # Compare complete direct truth-table transform, including vanished coefficients.
        n=2*m+t
        rows=[]
        for x in range(t):
            i,j=divmod(x,m)
            for y in range(1<<n):
                e,f,u=(y>>i)&1,(y>>(m+j))&1,(y>>(2*m+x))&1
                rows.append((-1)**(e^((e^f)&(e^u))))
        result=walsh(rows)
        decoded={beta:walsh(row) for beta,row in bands.items()}
        for index,v in enumerate(result):
            address,beta=divmod(index,1<<n)
            target=decoded.get(beta)
            assert v==(target[address]*(1<<(n-1)) if target else 0)
    return dict(m=m,blocks=t,s=s)

def unequal_clocks():
    p=[Fraction(x,8) for x in (4,2,1,1)]
    incident=[{0,2},{0,3},{1,2},{1,3}]
    appeared={(i,b):Fraction(0) for i,labels in enumerate(incident) for b in labels}
    for order in permutations(range(4)):
        weight=left=Fraction(1)
        seen=set()
        for i in order:
            weight*=p[i]/left
            left-=p[i]
        for i in order:
            for b in seen & incident[i]:
                appeared[i,b]+=weight
            seen.update(incident[i])
    for (i,b),prob in appeared.items():
        total=sum(p[j] for j in range(4) if b in incident[j])
        assert prob==1-p[i]/total
    return {f"{i},{b}":str(v) for (i,b),v in appeared.items()}

print(json.dumps(dict(python=platform.python_version(),grids=[grid(m) for m in (2,4,8)],
                      unequal_clock_probabilities=unequal_clocks()),indent=2))
