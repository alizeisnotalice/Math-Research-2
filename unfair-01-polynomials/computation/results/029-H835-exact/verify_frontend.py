"""Independent exact root and 32-bit necessary-layout front end for H8(3,5)."""
from collections import Counter
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import argparse
import gzip
import hashlib
import json
import struct

from exact_complex import cj, dis, lift, mother, pf


H = [1, 1, 1, 0, 0, 0, 0, 1, 1]
ALL = (1 << 32) - 1


def req(condition, message):
    if not condition:
        raise ValueError(message)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def local_energy(bits):
    return 5 * (bits & 7).bit_count() + 3 * (bits >> 3).bit_count()


def local_choices():
    choices = {}
    for bits in range(256):
        choices.setdefault(local_energy(bits), []).append(bits)
    return choices


def canonical_bytes():
    """Independently regenerate all quotient representatives in positional order."""
    choices = local_choices()
    req(choices[0] == [0] and choices[15] == [7, 248] and choices[30] == [255],
        'local terminal/trivial choices')
    req(all(local_energy(bits ^ 255) == 30 - local_energy(bits) for bits in range(256)),
        'local complement energy')
    out = bytearray()
    counts = Counter()
    for energy in range(1, 15):
        for blocks in product(choices.get(energy, ()), repeat=4):
            word = sum(block << (8 * index) for index, block in enumerate(blocks))
            out.extend(struct.pack('<I', word))
            counts[energy] += 1
    for tail in product(choices[15], repeat=3):
        word = 7 | tail[0] << 8 | tail[1] << 16 | tail[2] << 24
        out.extend(struct.pack('<I', word))
        counts[15] += 1
    expected = {3: 625, 5: 81, 6: 10000, 8: 50625, 9: 10000, 10: 81,
                11: 810000, 12: 625, 13: 50625, 14: 810000, 15: 8}
    req(counts == expected and sum(counts.values()) == len(out) // 4 == 1742670,
        'independent canonical counts')
    return bytes(out), counts


def disk(item):
    centre, radius = cj(item['disk']['center']), pf(item['disk']['radius'])
    req(radius > 0, 'positive disk radius')
    return centre, radius


def separated(disks):
    for index, (centre, radius) in enumerate(disks):
        for other, other_radius in disks[:index]:
            req(dis(centre, radius, other, other_radius), 'root disks disjoint')


def quadratic_box(centre, radius):
    def square_range(lower, upper):
        return (F(0) if lower <= 0 <= upper else min(lower * lower, upper * upper),
                max(lower * lower, upper * upper))
    real_square = square_range(centre[0] - radius, centre[0] + radius)
    imaginary_square = square_range(centre[1] - radius, centre[1] + radius)
    return [(real_square[0] + imaginary_square[0], real_square[1] + imaginary_square[1]),
            (-2 * (centre[0] + radius), -2 * (centre[0] - radius)), (F(1), F(1))]


def target_binary():
    left = [3 * index for index, coefficient in enumerate(H) if coefficient]
    right = [5 * index for index, coefficient in enumerate(H) if coefficient]
    coefficient = Counter(i + j for i in left for j in right)
    req(len(coefficient) == 25 and set(coefficient.values()) == {1} and max(coefficient) == 64,
        'fixed H8 target must be binary')


def verify(balls, ledger_path=None):
    req(balls['format'] == 'H8-a3b5-rational-rouche-v1' and balls['H'] == H and
        balls['a'] == 3 and balls['b'] == 5, 'fixed H8(3,5) target')
    target_binary()
    scale = int(balls['scale'])
    req(0 < scale <= 10 ** 1000, 'valid positive exact scale')
    req(len(balls['mothers']) == 4 and len(balls['groups']) == 32, 'H8 source cardinalities')

    mothers, mother_disks = {}, []
    for item in balls['mothers']:
        block = item['block']
        req(type(block) is int and 0 <= block < 4 and block not in mothers, 'unique mother block')
        centre, radius = disk(item)
        req(centre[1] > radius, 'positive nonreal mother representative')
        left, right = mother(H, centre, radius, scale)
        req(left < right and pf(item['lhs']) == left and pf(item['rhs']) == right, 'mother Rouché bound')
        mothers[block] = (centre, radius)
        mother_disks.extend(((centre, radius), ((centre[0], -centre[1]), radius)))
    req(len(mothers) == 4, 'all mother blocks')
    separated(mother_disks)

    seen, factors, target_disks = set(), {}, []
    family_degrees = Counter()
    for group in balls['groups']:
        index = group['index']
        req(type(index) is int and 0 <= index < 32 and index not in seen, 'unique group index')
        seen.add(index)
        expected_block, offset = divmod(index, 8)
        family, exponent, slot = ('a', 3, offset) if offset < 3 else ('b', 5, offset - 3)
        req(type(group['block']) is int and type(group['exponent']) is int and type(group['slot']) is int and
            group['block'] == expected_block and group['family'] == family and
            group['exponent'] == exponent and group['slot'] == slot and group['kind'] == 'quadratic',
            'group block/family/exponent/slot/kind binding')
        centre, radius = disk(group)
        req(abs(centre[1]) > radius, 'nonreal quadratic lift')
        left, right = lift(*mothers[expected_block], centre, radius, scale, exponent)
        req(left < right and pf(group['lhs']) == left and pf(group['rhs']) == right, 'lift Rouché bound')
        target_disks.extend(((centre, radius), ((centre[0], -centre[1]), radius)))
        factors[index] = quadratic_box(centre, radius)
        family_degrees[family] += 2
    req(len(seen) == len(factors) == 32 and family_degrees == Counter({'a': 24, 'b': 40}),
        'all H8 lift groups')
    req(len(target_disks) == 64, 'target roots through conjugates')
    separated(target_disks)
    req(balls['coverage'] == {'mother_pairs': 4, 'quadratic_groups': 32, 'root_disks': 64,
                               'pairwise_disjoint': True}, 'root coverage metadata')

    canonical, counts = canonical_bytes()
    if ledger_path is not None:
        ledger_path = Path(ledger_path)
        req(ledger_path.read_bytes() == canonical, 'external canonical u32le ledger exact ordering')
        packed = ledger_path.with_suffix(ledger_path.suffix + '.gz')
        if packed.is_file():
            req(gzip.decompress(packed.read_bytes()) == canonical, 'external gzip canonical ledger')
    standard = sum(7 << (8 * block) for block in range(4))
    req(sum(word[0] == standard for word in struct.iter_unpack('<I', canonical)) == 1,
        'unique aligned standard all-a word')
    result = {
        'status': 'independent_H835_exact_root_and_layout_frontend_passed',
        'degree': 64, 'binary_terms': 25, 'mother_roots': 8, 'target_roots': 64,
        'complex_parent_pairs': 4, 'quadratic_groups': 32, 'linear_groups': 0,
        'canonical_layouts': len(canonical) // 4,
        'counts_by_E': {str(key): counts[key] for key in sorted(counts)},
        'canonical_sha256': hashlib.sha256(canonical).hexdigest(),
        'standard_mask': standard,
        'scope': 'Necessary root-layout coverage only; no coefficient exclusions.',
    }
    return result, factors, canonical


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--balls', type=Path, required=True)
    parser.add_argument('--ledger', type=Path)
    parser.add_argument('--record', type=Path, required=True)
    args = parser.parse_args()
    result, _, _ = verify(json.loads(args.balls.read_text()), args.ledger)
    result.update(balls_sha256=sha(args.balls), verifier_sha256=sha(__file__),
                  helper_sha256=sha(Path(__file__).with_name('exact_complex.py')))
    if args.ledger:
        result['ledger_sha256'] = sha(args.ledger)
    args.record.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result))


if __name__ == '__main__':
    main()
