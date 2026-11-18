"""Construct exact closed boxes, including a tube connecting two verified low six-body poses."""
import argparse,json
from fractions import Fraction as F
from pathlib import Path

def make(c,root,kind):
    if kind=='tube':
        h0=list(map(F,json.loads((c/'free-p30-reoptimized-pose.json').read_bytes())[:11]))
        h1=[F(z,2**44) for z in json.loads((root/'other-pose.json').read_bytes())['center']]
        rows=[];N=192
        for j in range(N):
            row=[]
            for i,(a,b) in enumerate(zip(h0,h1)):
                lo,hi=sorted((a+(b-a)*F(j,N),a+(b-a)*F(j+1,N)))
                extra=F('.00002') if i in (2,5,8) else F('.000004')
                row.append([((lo-extra)*2**44).__floor__(),((hi+extra)*2**44).__ceil__()])
            rows.append(row)
        return rows
    rows=json.loads((root/'probe-poses.json').read_bytes())
    if kind=='residual':
        # Centered, closed, positive-volume portions of five reported first open boxes.
        return [[[(a+b)//2-(b-a)//(2*scale),(a+b)//2+(b-a)//(2*scale)] for a,b in row]
                for scale,row in zip((4,8,4,4,4),rows[:5])]
    if kind=='anchors':
        import sys
        from flint import ctx
        from hinge_bank import q
        sys.path.insert(0,str(c.resolve()))
        import certify_five_body_contact_stationary as st
        ctx.prec=192
        _,actual,_,_=st.model(*(st.interval(*b) for b in st.BOX))
        assert all(q(F(a,2**44))<=v<=q(F(b,2**44)) for (a,b),v in zip(rows[5],actual))
        pose=json.loads((c/'free-p30-reoptimized-pose.json').read_bytes())
        low=[[(F(v)*2**44).__floor__(),(F(v)*2**44).__ceil__()] for v in pose[:11]]
        return rows[5:]+[low]
    raise ValueError(kind)

if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('comp',type=Path);p.add_argument('kind',choices=['tube','anchors','residual']);p.add_argument('output',type=Path)
    a=p.parse_args();a.output.write_text(json.dumps(make(a.comp,Path(__file__).resolve().parent,a.kind))+'\n')
