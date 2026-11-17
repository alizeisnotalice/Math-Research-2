"""Direct exact AB verification of the displayed S2 last-row contradiction."""
import argparse,hashlib,json,math
from pathlib import Path
def main():
    if not __debug__:raise RuntimeError('no -O')
    ap=argparse.ArgumentParser();ap.add_argument('--witness',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);args=ap.parse_args();z=json.loads(args.witness.read_text())
    assert z['row']==[4,1,12,0,0,0,0,0] and len(z['GQ_at66'])==1
    G=z['GQ_at66'][0]['G'];Q=z['GQ_at66'][0]['Q'];assert len(G)==len(Q)==67 and G[0]==Q[0]==1 and Q[1:]==[0]*66
    h=hashlib.sha256((json.dumps([G,Q],separators=(',',':'))+'\n').encode()).hexdigest();assert h==z['state66_sha256']=='bf6c6fdb4b7680c0958c0979dbe48d783550899f08162a2b89e0d2aa06aeda87'
    def mul(x,y):return [x[0]*y[0]-x[1]*y[1],x[0]*y[1]+x[1]*y[0]+4*x[1]*y[1]]
    def sign(p,q):
        R=2*p+4*q
        if q==0:return (R>0)-(R<0)
        if q>0 and R>=0:return 1
        if q<0 and R<=0:return -1
        delta=R*R-12*q*q;assert delta!=0
        return ((delta>0)-(delta<0))*(1 if R>0 else -1)
    def unit(p,q):return sign(p,q)>=0 and sign(1-p,-q)>=0
    def factors(g,q):
        A=[];B=[]
        for k in range(len(g)):
            ag=g[k-2] if k>=2 else 0;aq=q[k-2] if k>=2 else 0
            A.append([g[k]+(g[k-7] if k>=7 else 0),ag])
            B.append([q[k]+(q[k-7] if k>=7 else 0)+4*aq,-aq])
        return A,B
    def coefficient(A,B,k):
        terms=[mul(A[i],B[k-i]) for i in range(k+1)];return [sum(x[j] for x in terms) for j in range(2)]
    A,B=factors(G,Q)
    assert all(unit(*x) for x in A+B)
    assert all(coefficient(A,B,k) in ([0,0],[1,0]) for k in range(67)) and A[5]==B[5]==[0,0] and G[3]==Q[3]==0
    g=-3787078993083270494;AA,BB=factors(G+[g],Q+[0]);p,q=AA[67]
    assert [p,q]==[-3668403729770544885,982945816903335992] and sign(p,q)>0 and sign(1-p,-q)>0
    # Integer shifts other than this g cannot stay in the unit interval.
    assert not unit(p-1,q) and not unit(p+1,q)
    got=[]
    for v in (0,1):
        A,B=factors(G+[g],Q+[v]);C=coefficient(A,B,67);assert C==[2+v,0]
        got.append({'g67':g,'q67':v,'A67':A[67],'B67':B[67],'C67':C})
    assert got==z['all67_choices']
    L=p+2*q;lo=3*q*q-L*L;hi=(1-L)**2-3*q*q;assert L<0 and q>0 and lo>0 and hi>0
    out={'status':'passed_all_67_prefix_rows_and_exhaustive_last_row','state66_sha256':h,'lower_square_margin':lo,'upper_square_margin':hi,'L':L,'q':q,'last_C_values':[2,3],'witness_sha256':hashlib.sha256(args.witness.read_bytes()).hexdigest(),'verifier_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    args.out.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out))
if __name__=='__main__':main()
