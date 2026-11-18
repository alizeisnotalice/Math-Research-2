"""Recover an actual improving cell and verify every Fourier coefficient on its bands."""
import hashlib
import json
import random

def recover(seed, expected_hash):
    n=12
    m=n//2
    size=1<<n
    width=1<<m
    leaves=size//2
    pts=list(range(size))
    random.Random(seed).shuffle(pts)
    pairs=list(zip(pts[::2],pts[1::2]))
    assert hashlib.sha256(json.dumps(pairs).encode()).hexdigest()==expected_hash
    assert len({x^y for x,y in pairs})>=2
    label={x:i for i,edge in enumerate(pairs) for x in edge}
    for answer in range(size//width):
        labels=[label[(answer<<m)+y] for y in range(width)]
        if len(set(labels))==width:
            break
    else:
        raise AssertionError('no certificate on this coordinate cut')
    magnitudes=set()
    for tag in labels:
        for a in range(width):
            c=sum((-1)**((a&y).bit_count()%2) for y,t in enumerate(labels) if t==tag)
            magnitudes.add(abs(c))
    assert magnitudes=={1}
    assert (width+m+1)*leaves>(leaves+n+1)*width
    triple=None
    for ans in range(8):
        active={}
        for y in range(512):
            active.setdefault(label[(y<<3)|ans],[]).append(y)
        f=sum(len(xs)==2 for xs in active.values())
        t=len(active)-f
        s3=256*f+512*t
        D3=len(active)+10
        if t and D3**2*leaves**2>(leaves+n+1)**2*s3:
            actual=sum(sum((-1)**((a&y).bit_count()%2) for y in xs)!=0
                       for xs in active.values() for a in range(512))
            assert actual==s3
            triple=dict(query='low three address bits',answer=ans,full_edges=f,
                        singles=t,D_child=D3,s_child=s3)
            break
    assert triple is not None
    return dict(seed=seed,three_query=triple,query='upper six address bits',answer=answer,
                distinct_labels=width,s_child=width**2,D_child=width+m+1,
                s_parent=leaves**2,D_parent=leaves+n+1)

if __name__=='__main__':
    proofs=[]
    for machine in ('mini','macbook'):
        with open('mixed-'+machine+'.json') as f:
            record=json.load(f)
        proofs.append(recover(record['seed'],record['edge_sha256']))
    print(json.dumps(proofs,indent=2))
