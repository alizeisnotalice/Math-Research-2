"""Close explicit positive-twist exterior residuals at the zero-slice threshold."""
from argparse import ArgumentParser
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction as Q
from hashlib import sha256
import json
from pathlib import Path
import sys


def main():
    parser=ArgumentParser()
    parser.add_argument('source',type=Path)
    parser.add_argument('stage',type=Path)
    parser.add_argument('output',type=Path)
    parser.add_argument('--computation',type=Path,default=Path(__file__).parent)
    parser.add_argument('--workers',type=int,default=2)
    parser.add_argument('--depth',type=int,default=12)
    args=parser.parse_args()
    sys.path.insert(0,str(args.computation))
    from certify_three_body_weighted import init_worker,solve

    raw=args.stage.read_bytes()
    records=json.loads(raw)
    assert records['target']=='0.834780945776'
    lo,hi=records['angle_indices']
    assert 294<=lo<hi<=1310 and tuple(map(Q,[records['frontier'][0][2],records['frontier'][1][2]]))==(Q(lo,65536),Q(hi,65536))
    jobs=[(0,i,tuple(tuple(map(Q,a)) for a in box)) for i,box in enumerate(records['residual'])]
    counts=Counter()
    residual=[]
    least=None
    ledger=sha256()
    with ProcessPoolExecutor(max_workers=args.workers,initializer=init_worker,initargs=(args.source,records['target'],args.depth)) as pool:
        for job,row in zip(jobs,pool.map(solve,jobs,chunksize=1)):
            assert row['index']==job[1]
            counts.update(row['counts'])
            residual.extend(row['unresolved'])
            if row['min_bound'] is not None:
                v=Q(row['min_bound']);least=v if least is None else min(least,v)
            ledger.update(json.dumps(row,sort_keys=True).encode()+b'\n')
    leaves=sum(v for key,v in counts.items() if key.startswith('weights_'))+len(residual)
    assert counts['checked']==2*leaves-len(jobs)
    result={'target':records['target'],'angle_indices':records['angle_indices'],'depth_limit':args.depth,'stage_sha256':sha256(raw).hexdigest(),
            'script_sha256':sha256(Path(__file__).read_bytes()).hexdigest(),
            'jobs':len(jobs),'counts':dict(counts),'remaining':residual,
            'weighted_min':str(least),'ledger_sha256':ledger.hexdigest()}
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print({k:v for k,v in result.items() if k!='remaining'},'remaining',len(residual),flush=True)


if __name__=='__main__':
    main()
