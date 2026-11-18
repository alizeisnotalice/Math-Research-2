"""Reconstruct every branch in trace files after replaying the numerical proof kernel."""
import argparse
from collections import Counter
from hashlib import sha256
import json
from pathlib import Path
import struct

WEIGHTS = (100, 100, 16, 100, 100, 11, 100, 100, 9, 100, 100)

def audit(root, trace):
    stack = [(root[:11], root[11:], 0)]
    cursor = visits = cuts = splits = depth = 0
    terminals = Counter()
    while stack:
        lo, hi, d = stack.pop()
        depth = max(depth,d)
        changed, kind, axis = struct.unpack_from("<Bbb", trace, cursor)
        cursor += 3
        visits += 1
        if kind:
            assert not changed and kind in (1, 2, 3)
            terminals[str(kind)] += 1
            continue
        assert 0 <= axis < 11
        if changed:
            ends = struct.unpack_from("<22q",trace,cursor)
            cursor += 176
            low, high = list(ends[:11]), list(ends[11:])
            assert all(a <= b <= c <= d for a,b,c,d in zip(lo,low,high,hi))
        else:
            low,high = lo,hi
        old = [(b-a)*w for a,b,w in zip(lo,hi,WEIGHTS)]
        new = [(b-a)*w for a,b,w in zip(low,high,WEIGHTS)]
        # The production kernel uses floating point only for its split and
        # prioritization decisions. Check independently with integral widths.
        if max((a-b)*10 for a,b in zip(old,new)) > max(100,max(old)):
            stack.append((low,high,d+1))
            cuts += 1
            continue
        if new[axis]*24 < max(new) or high[axis]-low[axis] < 2:
            axis = new.index(max(new))
        assert high[axis]-low[axis] >= 2
        mid = (low[axis]+high[axis])//2
        a,b = low.copy(),high.copy()
        a[axis] = b[axis] = mid
        stack.extend([(a,high,d+1),(low,b,d+1)])
        splits += 1
    assert cursor == len(trace) and sum(terminals.values()) == splits+1
    return dict(vis=visits,splits=splits,contractions=cuts,closed=dict(terminals),md=depth,remaining=[])

def main():
    if not __debug__:
        raise SystemExit("Run with Python assertions enabled.")
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("work",type=Path)
    p.add_argument("--expect",type=Path,help="directory of the published JSON records")
    a=p.parse_args()
    for pos in (139716,151617,139730,152529):
        stem=f"part1205-pos{pos}"
        f=a.work/(stem+".json")
        row=json.loads(f.read_bytes())
        if a.expect:
            assert row == json.loads((a.expect/f.name).read_bytes())
        trace=(a.work/(stem+".bin")).read_bytes()
        assert sha256(trace).hexdigest()==row["trace_sha256"]
        assert row["part"]==1205 and row["pos"]==pos and row["threshold"]=="0.83740"
        out=audit(row["root"],trace)
        assert out == {k:row["result"][k] for k in out}
        print(pos, out["vis"], sum(out["closed"].values()))

if __name__=="__main__":
    main()
