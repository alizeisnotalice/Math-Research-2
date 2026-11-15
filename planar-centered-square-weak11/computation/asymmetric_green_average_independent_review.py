#!/usr/bin/env python3
"""Independent 256-bit complex-root replay plus exact rectangle partition audit."""
from fractions import Fraction as Q
from pathlib import Path
import datetime
import hashlib
import json
import os
import platform
import subprocess
import sys
import time
import flint
from flint import arb,acb,ctx
ctx.prec=256


def fs(x):
    return str(x.numerator)+'/'+str(x.denominator)


def A(x):
    return arb(x.numerator)/x.denominator


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def endpoint_fraction(x,which):
    y=x.lower() if which=='lower' else x.upper()
    return Q(str(y.fmpq()))


def ball_record(x):
    lo=endpoint_fraction(x,'lower'); hi=endpoint_fraction(x,'upper')
    assert lo<=hi
    return {'ball':str(x),'lower_rational':fs(lo),'upper_rational':fs(hi)}


def inverse_root_green_log(s,b,point):
    # Direct complex quadratic root (sqrt(1+4*b*z/s)-1)/(2*b),
    # not the author's nested-real modulus formula or original evaluator.
    u=1/s
    if point=='plus' and u>=1+b:
        return arb(0)
    if point=='minus' and u>=1-b:
        return arb(0)
    zeta={'plus':acb(A(u)), 'minus':acb(-A(u)), 'imag':acb(0,A(u))}[point]
    if b==0:
        w=zeta
    else:
        w=((1+4*A(b)*zeta).sqrt()-1)/(2*A(b))
    value=-abs(w).log()
    assert value.is_finite()
    return value.max(arb(0))


def independent_bound(box):
    sl,sh,bl,bh=box
    lp=inverse_root_green_log(sh,bh,'plus')
    lm=inverse_root_green_log(sh,bl,'minus')
    li=inverse_root_green_log(sh,bh,'imag')
    upper_sum=lp+lm+2*li
    assert upper_sum>0
    numerator=arb.pi()*A(sl)**2*(1+A(bl)**2)
    bound=numerator/upper_sum
    return bound, {'L_plus':ball_record(lp),'L_minus':ball_record(lm),
                   'L_imag':ball_record(li),'denominator':ball_record(upper_sum)}


def partition_audit(boxes,root):
    rsl,rsh,rbl,rbh=root
    assert rsl<rsh and rbl<rbh
    for sl,sh,bl,bh in boxes:
        assert rsl<=sl<sh<=rsh and rbl<=bl<bh<=rbh
    overlaps=0
    for i,a in enumerate(boxes):
        for b in boxes[i+1:]:
            hit=max(a[0],b[0])<min(a[1],b[1]) and max(a[2],b[2])<min(a[3],b[3])
            overlaps+=hit
    assert overlaps==0
    s_cuts=sorted({rsl,rsh}|{x for box in boxes for x in box[:2]})
    b_cuts=sorted({rbl,rbh}|{x for box in boxes for x in box[2:]})
    count=0
    for sl,sh in zip(s_cuts,s_cuts[1:]):
        for bl,bh in zip(b_cuts,b_cuts[1:]):
            sm=(sl+sh)/2; bm=(bl+bh)/2
            multiplicity=sum(a<sm<c and d<bm<e for a,c,d,e in boxes)
            assert multiplicity==1,(fs(sm),fs(bm),multiplicity)
            count+=1
    total_area=sum((sh-sl)*(bh-bl) for sl,sh,bl,bh in boxes)
    root_area=(rsh-rsl)*(rbh-rbl)
    assert total_area==root_area
    return {'leaf_count':len(boxes),'pairwise_interior_overlap_count':overlaps,
            'pairs_checked':len(boxes)*(len(boxes)-1)//2,
            's_cut_count':len(s_cuts),'b_cut_count':len(b_cuts),
            'elementary_open_cells_with_multiplicity_one':count,
            'sum_leaf_area':fs(total_area),'root_area':fs(root_area),
            'closed_cover_reason':'Every open arrangement cell is covered once; finite union of closed rectangles covers the closure, including all shared and outer boundary lines. Only interiors are disjoint.'}


def main():
    assert os.path.ismount('/Volumes/data'), 'data volume not mounted'
    t0=time.monotonic()
    base=Path(__file__).resolve().parent.parent
    inp=base/'computation/inputs/asymmetric-green-average/vertex-leaves.json'
    data=json.loads(inp.read_text())
    root=tuple(map(Q,data['root_box']))
    boxes=[tuple(map(Q,row)) for row in data['leaf_boxes']]
    assert len(boxes)==102
    partition=partition_audit(boxes,root)
    threshold=Q(18,5)
    leaves=[]
    lowest=None
    for index,box in enumerate(boxes):
        bound,terms=independent_bound(box)
        record=ball_record(bound)
        lower=Q(record['lower_rational'])
        assert lower>threshold,(index,record)
        if lowest is None or lower<lowest[0]:
            lowest=(lower,index)
        leaves.append({'index':index,'box':[fs(x) for x in box],
                       'recomputed_necessary_bound':record,'terms':terms,
                       'strictly_above_18_over_5':True})
    tails={
        'small_s':9*arb.pi()/(64*(arb(9)/8).log()),
        'large_s':49*arb.pi()/(16*(arb(21)/4).log())}
    # These comparisons certify the monotone branches used in the analytic tails.
    assert (arb(9)/8).log()<arb(1)/2
    assert (arb(21)/4).log()>arb(1)/2
    for name,value in tails.items():
        assert endpoint_fraction(value,'lower')>threshold,(name,str(value))
    gitbase=subprocess.run(['git','-C',str(base.parent),'rev-parse','HEAD'],
                           capture_output=True,text=True,check=True).stdout.strip()
    review=base/'reviews/asymmetric-quadratic-green-average-independent-review.md'
    out={'status':'PASS_independent_102_leaf_256bit_complex_inverse_and_exact_partition',
         'timestamp_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
         'python':sys.version,'executable':sys.executable,'flint':flint.__version__,
         'platform':platform.platform(),'precision_bits':ctx.prec,'workers':1,
         'git_base':gitbase,'threshold':fs(threshold),
         'script_sha256':sha(Path(__file__).resolve()),'input_sha256':sha(inp),
         'review_sha256':sha(review) if review.exists() else None,
         'author_proof_sha256':data['author_proof_sha256'],
         'author_certificate_sha256':data['author_certificate_sha256'],
         'author_script_sha256':data['author_script_sha256'],
         'author_float_enclosures_used':False,'original_evaluator_imported':False,
         'expression':'Direct acb complex root; Green logs from -log(abs(inverse root)); all 102 boxes recomputed at 256 bits.',
         'partition':partition,'tail_bounds':{k:ball_record(v) for k,v in tails.items()},
         'minimum_central_leaf':{'index':lowest[1],'lower_rational':fs(lowest[0]),
                                 'lower_ball':leaves[lowest[1]]['recomputed_necessary_bound']['ball']},
         'leaves':leaves,'elapsed_seconds':time.monotonic()-t0,
         'scope':'Leaf coordinates only are reused as certificate input. Every cost enclosure is recomputed. The coverage and no-positive-area-overlap test uses exact Fraction arithmetic. Analytic geometry and monotonicity are reviewed separately. No claim about all non-D4 domains, R38 D4-domain mixtures, or actual full-boundary majorization.'}
    dest=base/'computation/results/asymmetric-green-average-independent-review/independent-256bit.json'
    dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:out[k] for k in ['status','python','flint','precision_bits','workers','partition','minimum_central_leaf','elapsed_seconds']},indent=2))

if __name__=='__main__':
    main()
