"""Check all 17 candidates using only signed integer evaluation."""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path

POINTS=((1,1),(1,2),(2,1),(1,4),(4,1),(3,4),(4,3),(6,5),(5,6))


def homogeneous(row,p,q,embedding):
    D,T_a,r,T_b,s,c,d,e,f,g=row
    tau=1 if D%4==1 else 2
    rational=(2,T_a,T_b,2*c,2*d,2*e,2*f,2*g,2)
    weights=[p**i*q**(8-i) for i in range(9)]
    rational_sum=sum(a*w for a,w in zip(rational,weights))
    return rational_sum, embedding*tau*(r*weights[1]+s*weights[2])


def sign(rational,radical,D):
    """Sign of rational+radical*sqrt(D), with a strict squared comparison."""
    if radical==0:
        return (rational>0)-(rational<0)
    if rational==0 or (rational>0)==(radical>0):
        return (radical>0)-(radical<0)
    square=rational*rational-D*radical*radical
    return ((square>0)-(square<0))*(1 if rational>0 else -1)


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('rows',type=Path)
    ap.add_argument('--out',type=Path,required=True)
    a=ap.parse_args()
    if not __debug__:
        raise RuntimeError('run without -O')
    payload=a.rows.read_bytes()
    rows=[json.loads(line) for line in payload.splitlines()]
    assert len(rows)==len({tuple(x[0]) for x in rows})==17
    ct=Counter()
    witnesses=[]
    for row,first_code in rows:
        D=row[0]
        evaluations=[(p,q,emb,homogeneous(row,p,q,emb)) for p,q in POINTS for emb in (1,-1)]
        failed=next((x for x in evaluations if sign(*x[3],D)<=0),None)
        if failed is not None:
            assert first_code=='root_positive_reject'
            p,q,embedding,value=failed
            ct[f'positive_{p}_{q}']+=1
        else:
            assert first_code=='root_outer_reject'
            p,q,embedding=-2,1,1
            value=homogeneous(row,p,q,embedding)
            # P(-2)<0 and real monic degree eight imply a root to the left of -2.
            assert sign(*value,D)<0
            ct['negative_axis_at_minus_two']+=1
        witnesses.append({'row':row,'p':p,'q':q,'embedding':embedding,'scaled_2q8_value':value})
    assert ct=={'positive_1_1':13,'positive_2_1':2,'negative_axis_at_minus_two':2}
    report={'status':'verified_all_direct_root_witnesses','binding_sha256':hashlib.sha256(payload).hexdigest(),
            'counts':dict(ct),'witnesses':witnesses}
    with a.out.open('x') as f:f.write(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report['counts']))


if __name__=='__main__':
    main()
