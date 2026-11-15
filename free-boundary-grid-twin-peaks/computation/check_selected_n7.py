"""Python-int audit of selected n=7 joint-peak counts.

The recurrence is evaluated only on subsets of each selected independent set;
the largest has 16 vertices, so this is independent from the full D4 table.
"""
from functools import cache
from math import factorial
from pathlib import Path
import json

N = 49

def neighbours(mask):
    result = 0
    for v in range(N):
        if not (mask >> v) & 1:
            continue
        x, y = v % 7, v // 7
        for X, Y in ((x-1,y), (x+1,y), (x,y-1), (x,y+1)):
            if 0 <= X < 7 and 0 <= Y < 7:
                result |= 1 << (Y*7 + X)
    return result

def joint_count(target):
    @cache
    def c(mask):
        if not mask:
            return factorial(N)
        numerator = sum(c(mask ^ (1 << v)) for v in range(N) if (mask >> v) & 1)
        denominator = mask.bit_count() + neighbours(mask).bit_count()
        q, r = divmod(numerator, denominator)
        assert r == 0
        return q
    return c(target)

def main():
    even_even = sum(1 << (7*y+x) for y in range(0,7,2) for x in range(0,7,2))
    targets = {
        'empty': 0,
        'corner': 1,
        'two_same_row': (1 << 0) | (1 << 2),
        'even_even_4x4': even_even,
    }
    result = json.loads((Path(__file__).parent / 'results/n7-exact-r005.json').read_text())
    actual = {name: joint_count(mask) for name, mask in targets.items()}
    assert {key: str(value) for key, value in actual.items()} == result['selected_joint_peak_counts']
    print(json.dumps({'verified': True, 'bit_lengths': {key: value.bit_length() for key, value in actual.items()},
                      'counts': {key: str(value) for key, value in actual.items()}}, indent=2))

if __name__ == '__main__':
    main()
