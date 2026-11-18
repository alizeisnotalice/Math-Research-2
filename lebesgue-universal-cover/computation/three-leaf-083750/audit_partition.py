"""Check the exact binary addresses joining strict exclusions and pending subtrees."""
from argparse import ArgumentParser
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
import sys

LEAF = '180/24225'
GOAL = F('.83750')


def read(path):
    return json.loads(path.read_bytes())


def bisect(bounds, steps):
    box = [list(r) for r in bounds]
    for code in steps:
        assert type(code) is int and 0 <= code < 22
        k, side = divmod(code, 2)
        lo, hi = box[k]
        assert lo < hi and (lo+hi)%2==0
        box[k][1-side] = (lo+hi)//2
    return box


def check_cover(root, receipts, frontiers, mode_class):
    used = set()
    unresolved = []
    counts = [0]*4
    nodes=0
    def walk(item):
        label = item['label']
        parent = item['steps']
        limits = item['bounds']
        assert item['leaf'] == LEAF and limits == bisect(root['bounds'], parent)
        if label not in receipts:
            unresolved.append(item)
            return
        nonlocal nodes
        row = receipts[label]
        assert label not in used
        used.add(label)
        assert row['leaf']==LEAF and F(row['goal'])==GOAL
        assert row['threshold_integer']==120696470013529293
        assert list(map(list,zip(row['lo'],row['hi'],strict=True)))==limits
        assert row['mode'] in mode_class
        stats = row['stats']
        assert len(stats)==35 and all(type(x) is int for x in stats)
        stops = stats[2:6]
        assert min(stops)>=0 and stats[1]>0
        if not row['mode']&1:
            assert stops[2]==0
        nodes+=stats[1]
        for i,s in enumerate(stops):
            counts[i]+=s
        p=frontiers.get((label,row['mode']),[])
        assert row['pending']==len(p)
        if not p:
            assert stats[0]==1 and stats[7]==0 and stats[1]==2*sum(stops)-1
            return
        assert stats[0]==0 and 1<=stats[7]<=3
        pending_visited=row.get('pending_visited',int(stats[7]!=1))
        assert 0<=pending_visited<=len(p)
        assert stats[1]==2*sum(stops)+len(p)-1+pending_visited
        assert list(map(list,zip(stats[11:22],stats[22:33],strict=True)))==p[0]['bounds']
        children=[x['steps'] for x in p]
        assert all(q[:len(parent)]==parent for q in children)
        assert all(len(q)>len(parent) for q in children) or (len(children)==1 and not sum(stops))
        sorted_paths=sorted(children)
        assert all(b[:len(c)]!=c for c,b in zip(sorted_paths,sorted_paths[1:]))
        axes={}
        for steps in children:
            for level in range(len(parent),len(steps)):
                prefix=tuple(steps[:level])
                axis=steps[level]//2
                assert axes.setdefault(prefix,axis)==axis
        for child in p:
            assert child['label'].startswith(label+'-')
            walk(child)
    walk(root)
    assert used==receipts.keys()
    residual_volume=sum((F(1,2**len(r['steps'])) for r in unresolved),F(0))
    assert 0<=residual_volume<=1
    return dict(complete=not unresolved, strict_threshold=str(GOAL), leaf=LEAF,
                nodes=nodes, terminals=dict(zip(('five','P30_cap','Q11_cap','mixed'),counts,strict=True)),
                stitched_runs=len(used), pending=len(unresolved), unresolved_volume_fraction=str(residual_volume),
                unresolved_volume_display=float(residual_volume), unresolved=unresolved)


def main(a):
    sys.path.insert(0,str(a.comp))
    from certify_tomar_six_high import FINE, SourceTree
    tree=SourceTree(a.tree)
    original=tree.leaf(*map(int,LEAF.split('/')))
    exact=[[x*FINE for x in r] for r in original]
    assert all(x.denominator==1 for r in exact for x in r)
    root=dict(leaf=LEAF,label=LEAF,bounds=[list(map(int,r)) for r in exact],steps=[])
    receipts={}
    assets={}
    for path in a.receipts:
        raw=path.read_bytes()
        assets[path.name]=sha256(raw).hexdigest()
        for line in raw.splitlines():
            item=json.loads(line)
            assert item['label'] not in receipts
            receipts[item['label']]=item
    frontiers={}
    for folder in a.frontier_dirs:
        for path in sorted(folder.glob('*-frontier.json')):
            raw=path.read_bytes()
            name=path.name
            label,mode=name[:-len('-frontier.json')].rsplit('-mode',1)
            label=label.replace('180-24225','180/24225',1)
            key=(label,int(mode))
            if key[0] in receipts and key[1]==receipts[key[0]]['mode']:
                data=json.loads(raw)
                assert key not in frontiers
                frontiers[key]=data
                assets[name]=sha256(raw).hexdigest()
    result=check_cover(root,receipts,frontiers,(0,2) if a.six_only else (0,1,2,3,4))
    result.update(assets=assets,original_bounds=root['bounds'])
    a.output.write_text(json.dumps(result,indent=2)+'\n')
    print({k:v for k,v in result.items() if k not in {'unresolved','assets','original_bounds'}},flush=True)


if __name__=='__main__':
    p=ArgumentParser(description=__doc__)
    for name in ('comp','tree','output'):
        p.add_argument(name,type=Path)
    p.add_argument('--receipts',nargs='+',type=Path,required=True)
    p.add_argument('--frontier-dirs',nargs='+',type=Path,required=True)
    p.add_argument('--six-only',action='store_true')
    main(p.parse_args())
