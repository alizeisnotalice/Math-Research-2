"""Count the legal local half-block cuts for every four-variable bent base."""
import itertools, multiprocessing, json
from collections import Counter
from pathlib import Path
def bent(mask):
    return all(sum((-1 if ((mask>>x)^(x&a).bit_count())&1 else 1) for x in range(16)) in (-4,4) for a in range(16))
def scan(args):
    bases,indices=args
    hist=Counter()
    for index in indices:
        legal={bases[index]^other for other in bases}
        for A in legal:
            if A.bit_count()!=8:continue
            points=[x for x in range(16) if A>>x&1]
            sums={x^y for x in points for y in points}
            cuts=0
            for subset in itertools.combinations(points,4):
                B=sum(1<<x for x in subset)
                if B in legal:
                    assert A^B in legal
                    cuts+=1
            hist['affine' if len(sums)==8 else 'nonaffine',cuts]+=1
    return hist
if __name__=='__main__':
    bases=[x for x in range(65536) if bent(x)]
    assert len(bases)==896
    with multiprocessing.Pool(12) as pool:
        parts=pool.map(scan,[(bases,range(i,len(bases),12)) for i in range(12)],1)
    hist=sum(parts,Counter());assert hist=={('affine',6):896*30,('nonaffine',2):896*360}
    out={str(k):v for k,v in hist.items()}
    Path(__file__).with_name('all_bases.json').write_text(json.dumps(out)+'\n');print(out)
