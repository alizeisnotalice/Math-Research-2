"""Check whole original leaves at 0.83744 with the scaled dense-material P30 guard."""
from argparse import ArgumentParser
from concurrent.futures import ThreadPoolExecutor
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
import subprocess
import sys

HERE = Path(__file__).resolve().parent


def digest(data):
    return sha256(data).hexdigest()


def main():
    if not __debug__:
        raise SystemExit('Assertions must be enabled.')
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('computation', type=Path)
    parser.add_argument('tree', type=Path)
    parser.add_argument('certificate', type=Path)
    parser.add_argument('--replay', nargs=3, type=Path, metavar=('TABLE64', 'WORK', 'CACHE'))
    parser.add_argument('--workers', type=int, default=2)
    args = parser.parse_args()
    assert 1 <= args.workers <= 32
    sys.path.insert(0, str(args.computation.resolve()))
    from certify_tomar_six_high import SourceTree, FINE

    cert = json.loads(args.certificate.read_bytes())
    assert Fraction(cert['target']) == Fraction('.83744')
    assert Fraction(cert['scale']) == Fraction('.999999006')
    goal = int(Fraction(cert['target']) * 2**57) + 1
    assert goal == 120687823102244742
    for name, expected_hash in cert['inputs'].items():
        assert digest((HERE / name).read_bytes()) == expected_hash
    tree = SourceTree(args.tree)
    keys = []
    visits = p30 = 0
    assert len(cert['receipts']) == len(cert['selections'])
    pairs = list(zip(cert['receipts'], cert['selections'], strict=True))
    for receipt, source in pairs:
        part, pos = receipt['part'], receipt['pos']
        keys.append((part, pos))
        assert receipt['complete'] and receipt['all_closed']
        target = Fraction(source['target'])
        assert Fraction(cert['target']) <= target < Fraction('.84')
        assert target * 2**57 != int(target * 2**57)
        assert receipt['target'] == str(target) and receipt['threshold_integer'] == int(target*2**57)+1
        expected = dict(part=part, pos=pos, target=source['target'], paths=[''], terminal_guard='mesh64-dense')
        assert source['driver'] in cert['inputs'] and source['driver'] in ('replay_scaled_p30_leaf_083744.py', 'replay_scaled_p30_high_leaf.py')
        assert source['driver'] != 'replay_scaled_p30_leaf_083744.py' or target == Fraction('.83744')
        original_bytes = source['selection_text'].encode('utf-8')
        assert json.loads(original_bytes) == expected and digest(original_bytes) == receipt['selection_sha256']
        box = tree.leaf(part, pos)
        lo, hi = [[coordinate[j] * FINE for coordinate in box] for j in (0, 1)]
        assert all(x.denominator == 1 for side in (lo, hi) for x in side)
        assert lo == receipt['original_lo'] and hi == receipt['original_hi']
        row, = receipt['rows']
        assert lo == row['lo'] and hi == row['hi']
        assert row['bits'] == '' and row['steps'] == [] and row['terminal_guard'] == 'mesh64-dense'
        st = row['stat']
        assert len(st) == 30 and st[0] == 1 and st[5] == 0 and st[1] == 2*(st[2]+st[3])-1
        assert st[1] >= 1 and st[2] >= 0 and st[3] >= 0 and 0 <= st[4] <= 240
        assert st[6] >= st[3] and st[29] <= st[3]
        visits += st[1]
        p30 += st[3]
    assert len(keys) == len(set(keys)) and keys == sorted(keys)
    assert list(map(list, keys)) == cert['original_leaves']
    assert (visits, p30) == (cert['visits'], cert['p30_terminals'])
    if args.replay:
        table, work, cache = args.replay
        work.mkdir(parents=True, exist_ok=True)
        def run(pair):
            receipt, source = pair
            p, i = receipt['part'], receipt['pos']
            group = work / f'{p}-{i}-selection.json'
            group.write_bytes(source['selection_text'].encode('utf-8'))
            dest = work / f'{p}-{i}'
            subprocess.run([sys.executable, '-I', '-B', HERE / source['driver'],
                            args.computation, args.tree, table, dest, cache, group], check=True)
            assert json.loads((dest / 'result.json').read_bytes()) == receipt
        with ThreadPoolExecutor(args.workers) as pool:
            list(pool.map(run, pairs))
    print(json.dumps(dict(complete=True, material_scale=cert['scale'], target=cert['target'],
                          leaves=len(keys), visits=visits, p30_terminals=p30,
                          geometric_replayed=bool(args.replay))))


if __name__ == '__main__':
    main()
