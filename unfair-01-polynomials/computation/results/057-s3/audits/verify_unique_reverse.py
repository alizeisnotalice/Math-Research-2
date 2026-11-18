"""Expand the first empty layer for the unique S3 mother needing reversal."""
import argparse
import json
from pathlib import Path
import s3_trees as exact


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--out',type=Path,required=True)
    a=p.parse_args()
    mother=(0,0,7,1,45,1,0,0,0)
    P,Pb,T,N,D=exact.setup(mother,True)
    states=[((1,),(1,),(exact.ONE,),(exact.ONE,))]
    records=[]
    for k in range(1,27):
        next_states=[]
        attempts=[]
        for G,Q,A,B in states:
            base_a=exact.cofactor_base(P,G,k)
            base_b=exact.cofactor_base(Pb,Q,k)
            cross=exact.ZERO
            for i in range(1,k):
                cross=exact.plus(cross,exact.times(A[i],B[k-i],T,N))
            left=exact.shifts(base_a,T,D,k,True)
            right=exact.shifts(base_b,T,D,k,True)
            for g,ak in left:
                for q,bk in right:
                    val=exact.plus(cross,exact.plus(ak,bk))
                    if k==26:
                        attempts.append({'G':G,'Q':Q,'base_a':base_a,'base_b':base_b,'g26':g,'q26':q,'C26':val})
                    if val in (exact.ZERO,exact.ONE):
                        next_states.append((G+(g,),Q+(q,),A+(ak,),B+(bk,)))
        records.append(len(next_states))
        states=next_states
    assert records[-1]==0 and min(records[:-1])>0 and max(records)==237 and attempts
    result={'status':'verified_unique_reverse_empty_layer','mother':mother,'counts':records,'last_layer_attempts':attempts}
    with a.out.open('x') as f:f.write(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'status':result['status'],'parents_at_layer_25':records[-2],'attempts_at_layer_26':len(attempts),'values':sorted({tuple(a['C26']) for a in attempts})}))


if __name__=='__main__':
    main()
