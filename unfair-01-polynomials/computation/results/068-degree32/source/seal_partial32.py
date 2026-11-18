"""Preserve flushed proof rows from a generation shard stopped for rebalancing."""
from argparse import ArgumentParser
from collections import Counter
import json
from pathlib import Path

from census32 import sha


def encode(item):
    return json.dumps(item,sort_keys=True,separators=(',',':')).encode()+b'\n'


def main():
    if not __debug__:
        raise SystemExit('assertions required')
    p=ArgumentParser(description=__doc__)
    p.add_argument('--census',type=Path,required=True)
    p.add_argument('--partial',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--receipt',type=Path,required=True)
    p.add_argument('--exit-code',type=int,choices=(143,),required=True)
    a=p.parse_args()
    counts=Counter()
    candidates=json.loads(a.census.read_bytes())['ledger']['candidates']
    with a.partial.open('rb') as src:
        header=json.loads(next(src))
        assert header['mode']=='generate' and header['census_sha256']==sha(a.census)
        lo,hi=header['start'],header['stop']
        n=lo
        lines=[]
        discarded_bytes=0
        for raw in src:
            if not raw.endswith(b'\n'):
                discarded_bytes=len(raw)
                assert not src.read(1)
                break
            row=json.loads(raw)
            assert row.get('index')==n and row['shape']==candidates[n]['shape']
            counts[row['status']]+=1
            if row['status']=='unit':
                counts['mod'+str(row['prime'])]+=1
                counts[row['method']]+=1
            else:
                assert row['status']=='unresolved'
            lines.append(raw)
            n+=1
    assert lo<n<hi and counts['unit']+counts['unresolved']==n-lo
    interruption={'raw_partial_sha256':sha(a.partial),'generator_exit_code':a.exit_code,
                  'start':lo,'planned_stop':hi,'flushed_stop':n,'discarded_bytes':discarded_bytes,
                  'source_sha256':sha(Path(__file__))}
    header.update({'stop':n,'checkpoint':interruption})
    tail={'status':'complete','start':lo,'stop':n,'seconds':None,'counts':dict(counts)}
    with a.out.open('xb') as f:
        f.write(encode(header))
        f.writelines(lines)
        f.write(encode(tail))
    interruption['reframed_sha256']=sha(a.out)
    with a.receipt.open('xb') as f:
        f.write(encode(interruption))
    print(n,flush=True)


if __name__=='__main__':
    main()
