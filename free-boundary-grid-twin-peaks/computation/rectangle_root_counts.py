"""Python arbitrary-int independent-set cross-check for rectangle root bands."""
from math import factorial, comb
from fractions import Fraction
import json


def rectangle(rows,cols):
    N=rows*cols;valid=[s for s in range(1<<cols) if not s&(s<<1)];bags=[[] for _ in range(N+1)]
    def generate(row,previous,mask):
        if row==rows:bags[mask.bit_count()].append(mask);return
        for s in valid:
            if not s&previous:generate(row+1,s,mask|(s<<(row*cols)))
    generate(0,0,0)
    closed=[]
    for v in range(N):
        x,y=v%cols,v//cols;m=1<<v
        for X,Y in ((x-1,y),(x+1,y),(x,y-1),(x,y+1)):
            if 0<=X<cols and 0<=Y<rows:m|=1<<(Y*cols+X)
        closed.append(m)
    band=sum(1<<(y*cols+x) for y in range(rows) for x in range((cols-1)//2+1))
    C={0:factorial(N)};closures={0:0};L=0;A1=0
    for k in range(1,N+1):
        for s in bags[k]:
            b=s&-s;cl=closures[s^b]|closed[b.bit_length()-1];closures[s]=cl
            total=0;q=s
            while q:
                b=q&-q;q-=b;total+=C[s^b]
            c,rem=divmod(total,cl.bit_count());assert rem==0;C[s]=c
            sign=1 if k%2 else -1;L+=sign*(s&band).bit_count()*c;A1+=sign*k*c
    assert 0<L<=A1<=factorial(N)
    if cols%2==0:assert 2*L==A1
    return {'rows':rows,'cols':cols,'independent_sets':sum(map(len,bags)),'A1':A1,'L':L}

rows=[rectangle(*p) for p in [(5,2),(5,3),(6,2),(6,3),(6,4)]]
expected=[38272,2409332016,973824,1409237297920,2669227667409855360]
assert [z['L'] for z in rows]==expected
L={z['cols']:z['L'] for z in rows if z['rows']==6};A2=304351937651819767235871772704107520
C={j:comb(36,6*j)*L[j]*L[6-j] for j in (2,3,4)}
print(json.dumps({'verified':True,'rectangles':rows,'n6_cut_counts':C,
 'n6_half_sum_lower_bound':str(Fraction(sum(C.values()),2*A2)),
 'n6_parity_lower_bound':str(Fraction(max(C[2]+C[4],C[3]),A2))},indent=2))
