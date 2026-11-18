"""Rebuild the full degree-31 census through the older tree and tail pilot."""
from argparse import ArgumentParser
from collections import Counter
from hashlib import file_digest, sha256
from itertools import product
import json
from pathlib import Path
import sys
import time

sys.dont_write_bytecode=True
from census31 import OLD, SOURCE_NAMES, sha
from audit_degree31 import tail_prefix
from prefix_29_30 import early_code_tree


def main():
    if not __debug__:
        raise SystemExit('assertions required')
    parser=ArgumentParser(description=__doc__)
    parser.add_argument('--census',type=Path,required=True)
    parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args()
    stamp=time.monotonic()
    original=json.loads(args.census.read_bytes())
    seen={row['shape']:(status,row) for status,rows in original['ledger'].items() for row in rows}
    assert len(seen)==27648
    counts=Counter()
    normalized=sha256()
    for row in product('ABX','ABX','ABX','ABXZ','ABXZ','ABXZ','ABXZ','ABXZ'):
        word=''.join(row)
        prior,entry=seen[word]
        terms=list(early_code_tree(31,word,None,[(1,)]*12,Counter()))
        if not terms:
            status='early_rejections'
        else:
            assert len(terms)==1
            shape,A,B=terms[0]
            assert shape==(word,None,(1,)*12)
            tail=tail_prefix(A,B,(1,)*7)
            status='late_rejections' if tail is None else 'candidates'
            if tail is not None:
                left,right=tail
                assert entry['A']==[[x.lo,x.hi] for x in left]
                assert entry['B']==[[x.lo,x.hi] for x in right]
        assert prior==status
        normalized.update(f'{word} {status}\n'.encode())
        counts[status]+=1
    assert dict(counts)=={k:v for k,v in original['counts'].items() if k!='all'}
    here=Path(__file__).resolve().parent
    source={k:sha(OLD/k) for k in SOURCE_NAMES}
    source.update({n:sha(here/n) for n in ('audit_census31.py','audit_degree31.py','census31.py')})
    record={'status':'verified','counts':dict(counts),'census_sha256':sha(args.census),
            'ordered_shape_status_sha256':normalized.hexdigest(),'sources_sha256':source,
            'python':sys.version,'seconds':time.monotonic()-stamp}
    with args.out.open('x') as out:
        out.write(json.dumps(record,sort_keys=True,indent=2)+'\n')
    print(record['counts'],record['seconds'],flush=True)


if __name__=='__main__':
    main()
