"""Independent necessary-mask reconstruction for H8(x^3)H8(x^5).

This is deliberately a layout audit only.  It never evaluates factor
coefficients and hence cannot exclude a layout or settle the general problem.
"""
from collections import Counter
from itertools import product
from pathlib import Path
import gzip
import hashlib
import json
import struct


ROOT = Path('/Volumes/data/CodexTemp/T2/unfair-01-polynomials')
INPUT = ROOT / 'agents/p8_next_input/H8_a3b5_layout'
LEDGER = INPUT / 'H8-a3b5-canonical-masks.u32le'
PACKED = INPUT / 'H8-a3b5-canonical-masks.u32le.gz'
SUMMARY = INPUT / 'H8-a3b5-layout-summary.json'
HERE = Path(__file__).resolve().parent

H = [1, 1, 1, 0, 0, 0, 0, 1, 1]
ALL = (1 << 32) - 1


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def local_energy(bits):
    """One conjugate parent pair: three a and five b quadratic lifts."""
    return 5 * (bits & 0b111).bit_count() + 3 * (bits >> 3).bit_count()


def local_choices():
    answer = {}
    for bits in range(256):
        answer.setdefault(local_energy(bits), []).append(bits)
    return answer


def pack(blocks):
    assert len(blocks) == 4
    return sum(value << (8 * index) for index, value in enumerate(blocks))


def canonical_words():
    """The complete quotient representatives, regenerated without ledger input."""
    choices = local_choices()
    # E<15 selects the lower member of each complement pair.
    for energy in range(1, 15):
        for blocks in product(choices.get(energy, ()), repeat=4):
            yield energy, pack(blocks)
    # At E=15 exactly one member of each complement pair has block zero all-a.
    all_a = 0b00000111
    assert local_energy(all_a) == 15
    for tail in product(choices[15], repeat=3):
        yield 15, pack((all_a,) + tail)


def word_check(energy, word):
    assert type(word) is int and 0 <= word <= ALL
    blocks = tuple((word >> (8 * index)) & 0xff for index in range(4))
    assert tuple(local_energy(block) for block in blocks) == (energy,) * 4
    complement = word ^ ALL
    complement_blocks = tuple((complement >> (8 * index)) & 0xff for index in range(4))
    assert tuple(local_energy(block) for block in complement_blocks) == (30 - energy,) * 4
    if energy == 15:
        assert blocks[0] == 0b00000111
        assert complement_blocks[0] == 0b11111000


def target_check():
    left = [3 * degree for degree, coefficient in enumerate(H) if coefficient]
    right = [5 * degree for degree, coefficient in enumerate(H) if coefficient]
    convolution = Counter(a + b for a in left for b in right)
    assert len(convolution) == 25
    assert set(convolution.values()) == {1}
    assert min(convolution) == 0 and max(convolution) == 64
    return {'left_support': left, 'right_support': right,
            'support': sorted(convolution), 'binary_terms': len(convolution)}


def main():
    target = target_check()
    data = LEDGER.read_bytes()
    assert len(data) % 4 == 0
    assert gzip.decompress(PACKED.read_bytes()) == data
    choices = local_choices()
    assert choices[15] == [0b00000111, 0b11111000]

    counts = Counter()
    cursor = 0
    standard = pack((7, 7, 7, 7))
    standard_seen = 0
    for position, (energy, rebuilt) in enumerate(canonical_words()):
        word = struct.unpack_from('<I', data, cursor)[0]
        cursor += 4
        assert word == rebuilt, f'ledger differs at canonical position {position}'
        word_check(energy, word)
        counts[energy] += 1
        if word == standard:
            assert energy == 15
            standard_seen += 1
    assert cursor == len(data)
    assert standard_seen == 1
    total = sum(counts.values())
    assert total == len(data) // 4

    # This is an algebraic enumeration count, rather than an imported count.
    derived_counts = {}
    for energy, values in choices.items():
        if 0 < energy < 15:
            derived_counts[energy] = len(values) ** 4
    derived_counts[15] = len(choices[15]) ** 3
    assert counts == Counter(derived_counts)

    source_summary = json.loads(SUMMARY.read_text())
    assert source_summary['target']['H_coefficients_ascending'] == H
    assert source_summary['groups']['quadratic_group_count'] == 32
    assert source_summary['groups']['linear_group_count'] == 0
    assert source_summary['ledger']['uncompressed_sha256'] == sha(LEDGER)
    assert source_summary['ledger']['compressed_sha256'] == sha(PACKED)

    receipt = {
        'status': 'independent_H835_32bit_canonical_layout_audit_passed',
        'target': target,
        'groups': {'complex_parent_pairs': 4, 'quadratic_groups': 32, 'linear_groups': 0,
                   'per_block': {'a': 3, 'b': 5}},
        'balance': 'E=5*u+3*v per block; common E across all four blocks',
        'canonical_rule': 'E<15, or E=15 with first block all-a (0b00000111)',
        'counts_by_E': {str(key): counts[key] for key in sorted(counts)},
        'canonical_mask_count': total,
        'standard_mask': standard,
        'standard_mask_count': standard_seen,
        'ledger_sha256': sha(LEDGER),
        'compressed_sha256': sha(PACKED),
        'source_summary_sha256': sha(SUMMARY),
        'auditor_sha256': sha(Path(__file__)),
        'scope': 'Exact necessary root-layout coverage only; no root disks or coefficient-sign scan has run.',
    }
    out = HERE / 'H835-layout-audit.json'
    out.write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps(receipt))


if __name__ == '__main__':
    main()
