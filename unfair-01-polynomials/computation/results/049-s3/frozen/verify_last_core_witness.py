"""Exact independent check of the unique depth-52 prefix and row-53 conflict."""
import argparse,hashlib,json
from pathlib import Path

def sign(x):
    p,q=x;L=2*p+5*q
    if not q:return (L>0)-(L<0)
    if L>=0 and q>0:return 1
    if L<=0 and q<0:return -1
    if L>0:return (L*L>17*q*q)-(L*L<17*q*q)
    return (17*q*q>L*L)-(17*q*q<L*L)
def mul(x,y):return x[0]*y[0]-2*x[1]*y[1],x[0]*y[1]+x[1]*y[0]+5*x[1]*y[1]
def coeff(p,g,k):return (sum(p[i][0]*g[k-i] for i in range(8) if 0<=k-i<len(g)),sum(p[i][1]*g[k-i] for i in range(8) if 0<=k-i<len(g)))
def convolution(A,B,k):
    terms=[mul(A[i],B[k-i]) for i in range(k+1)]
    return sum(x[0] for x in terms),sum(x[1] for x in terms)

def main():
    if not __debug__:raise RuntimeError('no -O')
    ap=argparse.ArgumentParser();ap.add_argument('--witness',type=Path,required=True);ap.add_argument('--tree',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args()
    w=json.loads(a.witness.read_text());tree=json.loads(a.tree.read_text());assert w['row']==tree['row']==[5,2,17,0,0,0,0]
    G,Q=w['G'],w['Q'];assert len(G)==len(Q)==53 and Q==[1]+[0]*52
    digest=hashlib.sha256(json.dumps([G,Q],separators=(',',':')).encode()+b'\n').hexdigest()
    assert tree['trace'][51]['k']==52 and tree['trace'][51]['live']==1 and tree['trace'][51]['state_sha256']==digest
    P=[(1,0),(0,0),(0,0),(0,1),(0,0),(0,0),(0,0),(1,0)];Pb=[(u+5*v,-v) for u,v in P]
    A=[coeff(P,G,k) for k in range(53)];B=[coeff(Pb,Q,k) for k in range(53)]
    for z in A+B:assert sign(z)>=0 and sign((z[0]-1,z[1]))<=0
    assert all(convolution(A,B,k) in ((0,0),(1,0)) for k in range(53))
    pg=coeff(P,G,53);pq=coeff(Pb,Q,53);assert list(pg)==w['A53_before_free_integer'] and pq==(0,0)
    records=w['all_integer_translations'];assert len(records)==2
    # For the irrational row the interval length is one with irrational
    # endpoints, so any admissible integer shift is unique. The rational row
    # is exactly q53 in {0,1}. Thus checking this list proves exhaustion.
    g=records[0]['G53'];z=(pg[0]+g,pg[1]);assert pg[1]!=0 and sign(z)>0 and sign((z[0]-1,z[1]))<0
    values=[]
    for q,rec in enumerate(records):
        assert rec['G53']==g and rec['Q53']==q
        aa=A+[z];bb=B+[(q,0)];value=convolution(aa,bb,53);assert value==(2+q,0) and list(value)==rec['C53'];values.append(value)
    assert tree['status']=='dead_all_prefixes' and tree['depth_reached']==53
    out={'status':'passed_complete_prefix_and_exhaustive_row53_conflict','C53':values,'scope':'Uniqueness at row52 is bound to the separately fully replayed tree; this checker independently verifies the displayed coefficients and exhaustive row53 choices.','witness_sha256':hashlib.sha256(a.witness.read_bytes()).hexdigest(),'tree_sha256':hashlib.sha256(a.tree.read_bytes()).hexdigest(),'verifier_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    a.out.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out))

if __name__=='__main__':main()
