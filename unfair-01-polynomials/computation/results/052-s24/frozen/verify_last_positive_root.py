"""Exact analytic exclusion of the sole persistent S24 prefix mother."""
import argparse,hashlib,json
from fractions import Fraction as F
from pathlib import Path
def main():
    if not __debug__:raise RuntimeError('no -O')
    ap=argparse.ArgumentParser();ap.add_argument('--certificate',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args();z=json.loads(a.certificate.read_text());assert z['row']==[13,1,1,-2,1,3,-1,0]
    x=F(6,5);coeff=[(F(1),F(0)),(F(0),F(0)),(F(3,2),F(-1,2)),(F(3),F(0)),(F(-3,2),F(-1,2)),(F(-1),F(0)),(F(0),F(0)),(F(1),F(0))]
    p=sum(u*x**i for i,(u,v) in enumerate(coeff));q=sum(v*x**i for i,(u,v) in enumerate(coeff));margin=13*q*q-p*p
    assert p==F(494411,78125) and q==F(-1098,625) and p>0 and q<0 and margin==F(446075579,6103515625)>0
    assert z['point']=='6/5' and z['Pbar_value']==[str(p),str(q)] and z['negative_sign_square_margin']==str(margin)
    out={'status':'passed_exact_Pbar_positive_root_exclusion','row':z['row'],'point':str(x),'p':str(p),'q':str(q),'strict_negative_square_margin':str(margin),'certificate_sha256':hashlib.sha256(a.certificate.read_bytes()).hexdigest(),'verifier_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'proof':'Pbar(0)=1, Pbar(6/5)<0 implies a root in (0,6/5). Pbar divides C. Every nonnegative C with C(0)=1 is positive for x>0, contradiction.'}
    a.out.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out))
if __name__=='__main__':main()
