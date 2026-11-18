"""Exact support and dyadic interval census for the fixed all-one rows through 30."""
from argparse import ArgumentParser
from collections import Counter
from hashlib import file_digest
from itertools import product
import json
from pathlib import Path
import sys
import time

sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
OLD=HERE.parent/'round065-gap7-next'
sys.path.insert(0,str(OLD))
from prefix_29_30 import fill_interval, point, t_box
from verify_degrees_17_18 import common_values, real_bounds
from verify_degrees_25_26 import setup

DEGREE=31
PARTNER={'A':'A','B':'B','X':'Z','Z':'F'}
TAIL_RULE='FFBAABZ'
EARLY_BITS=(1,)*12
TAIL_BITS=(1,)*7
SOURCE_NAMES=('prefix_29_30.py','verify_degree_19.py','verify_degrees_17_18.py',
              'verify_degrees_21_22.py','verify_degrees_23_24.py','verify_degrees_25_26.py')


def sha(path):
    with path.open('rb') as f:
        return file_digest(f,'sha256').hexdigest()


def descriptor(early):
    return early,None,EARLY_BITS


def shapes():
    for v in product('ABX','ABX','ABX',*('ABXZ',)*5):
        yield ''.join(v)


def interval_state(early):
    rules=early+''.join(PARTNER[code] for code in early[::-1])+TAIL_RULE
    assert len(rules)==23
    A,B=map(lambda side:[point(x) for x in side],setup(t_box()))
    for row,state in enumerate(rules,8):
        interior=sum(a*b for a,b in zip(A[1:],reversed(B[1:]),strict=True))
        remaining=1-interior
        strict=row<16 and state!='Z'
        pair=fill_interval(remaining,state,strict)
        if pair is None:
            return {'shape':early,'failure_row':row,'support':state,
                    'remaining':[remaining.lo,remaining.hi]}
        A.append(pair[0]); B.append(pair[1])
    return {'shape':early,'A':[[x.lo,x.hi] for x in A],'B':[[x.lo,x.hi] for x in B]}


def census():
    real_bounds()
    start=time.monotonic()
    report={'early_rejections':[],'late_rejections':[],'candidates':[]}
    failed_by_row=Counter()
    for word in shapes():
        row=interval_state(word)
        failed=row.get('failure_row')
        if failed is None:
            report['candidates'].append(row)
        else:
            key='early' if failed<24 else 'late'
            report[key+'_rejections'].append(row)
            failed_by_row[failed]+=1
    assert {k:len(v) for k,v in report.items()}=={'early_rejections':sum(n for k,n in failed_by_row.items() if k<24),
                                                  'late_rejections':sum(n for k,n in failed_by_row.items() if k>=24),
                                                  'candidates':len(report['candidates'])}
    counts={k:len(v) for k,v in report.items()}
    counts['all']=sum(counts.values())
    assert counts['all']==3**3*4**5==27648
    src={n:sha(OLD/n) for n in SOURCE_NAMES}
    src[Path(__file__).name]=sha(Path(__file__))
    return {'status':'complete','degree':DEGREE,'counts':counts,'failure_rows':dict(sorted(failed_by_row.items())),
            'early_word':EARLY_BITS,'late_word':TAIL_BITS,'ledger':report,
            'python':sys.version,'seconds':time.monotonic()-start,'sources_sha256':src}


def main():
    if not __debug__:
        raise SystemExit('assertions required')
    pa=ArgumentParser(description=__doc__)
    pa.add_argument('--out',type=Path,required=True)
    args=pa.parse_args()
    result=census()
    with args.out.open('x') as out:
        out.write(json.dumps(result,sort_keys=True,separators=(',',':'))+'\n')
    print(result['counts'],result['failure_rows'],result['seconds'],flush=True)


if __name__=='__main__':
    main()
