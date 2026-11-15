"""Independent finite certificate for an all-n row/column-orbit obstruction.

Checks all 5!^2 intact core permutations. Only the four K=2 cases need
single-cut checks: deleting edges can only add peaks.
"""
from itertools import permutations
import json
CORE=(25,7,24,4,16,9,13,23,20,22,21,6,12,15,19,3,1,14,17,5,10,8,11,2,18)

def peaks(a,mask=0):
    out=[]
    for v in range(25):
        x,y=v%5,v//5;nn=[]
        if x and not(mask>>(x-1)&1):nn.append(v-1)
        if x<4 and not(mask>>x&1):nn.append(v+1)
        if y and not(mask>>(4+y-1)&1):nn.append(v-5)
        if y<4 and not(mask>>(4+y)&1):nn.append(v+5)
        if all(a[v]>a[z] for z in nn):out.append(v)
    return out

ps=list(permutations(range(5)));minimum=25;two=[]
for rows in ps:
    for cols in ps:
        a=tuple(CORE[y*5+x] for y in rows for x in cols);pp=peaks(a);minimum=min(minimum,len(pp))
        assert len(pp)>=2
        if len(pp)==2:
            d=abs(pp[0]%5-pp[1]%5)+abs(pp[0]//5-pp[1]//5)
            assert d==2
            cut_counts=[len(peaks(a,1<<j)) for j in range(8)]
            assert min(cut_counts)>=3
            two.append({'rows':rows,'cols':cols,'labels':a,'peaks':pp,'single_cut_peak_counts':cut_counts})
assert len(two)==4 and minimum==2
# Choose an actual K=2 core ordering and monotonically attach the low exterior.
def padded(n):
    low=n*n-25;a=[0]*(n*n)
    for y in range(5):
        for x in range(5):a[y*n+x]=low+two[0]['labels'][y*5+x]
    outside=[(max(x-4,0)+max(y-4,0),y,x) for y in range(n) for x in range(n) if x>=5 or y>=5]
    for label,(_,y,x) in zip(range(low,0,-1),sorted(outside)):a[y*n+x]=label
    pp=[]
    for y in range(n):
        for x in range(n):
            ns=[(Y,X) for Y,X in ((y-1,x),(y+1,x),(y,x-1),(y,x+1)) if 0<=Y<n and 0<=X<n]
            if all(a[y*n+x]>a[Y*n+X] for Y,X in ns):pp.append(y*n+x)
    assert sorted(a)==list(range(1,n*n+1)) and len(pp)==2
    assert abs(pp[0]%n-pp[1]%n)+abs(pp[0]//n-pp[1]//n)==2
    return {'n':n,'peaks':pp,'distance':2}
print(json.dumps({'verified':True,'core':CORE,'intact_permutations':len(ps)**2,'minimum_peaks':minimum,
                  'two_peak_permutations':two,'padded_witnesses':[padded(n) for n in (5,6,12,32)]},indent=2))
