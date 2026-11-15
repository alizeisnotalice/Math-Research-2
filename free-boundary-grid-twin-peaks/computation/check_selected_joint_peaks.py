"""Independent Python-int checks of selected 6x6 joint-peak counts.

It recomputes the recurrence only on every subset of each selected independent
set.  Thus it has different integer arithmetic and never uses the C++ table.
"""
from functools import cache
import json
from math import factorial
from pathlib import Path
import subprocess
import sys


N = 36
NFACT = factorial(N)


def neighbours(mask: int) -> int:
    answer = 0
    for v in range(N):
        if mask >> v & 1:
            x, y = v % 6, v // 6
            for X, Y in ((x-1, y), (x+1, y), (x, y-1), (x, y+1)):
                if 0 <= X < 6 and 0 <= Y < 6:
                    answer |= 1 << (Y*6+X)
    return answer


def joint_count(target: int) -> int:
    @cache
    def c(mask: int) -> int:
        if not mask:
            return NFACT
        numerator = sum(c(mask ^ (1 << v)) for v in range(N) if mask >> v & 1)
        denominator = mask.bit_count() + neighbours(mask).bit_count()
        quotient, remainder = divmod(numerator, denominator)
        assert remainder == 0
        return quotient
    return c(target)


def main(binary: Path):
    even = sum(1 << (y*6+x) for y in range(6) for x in range(6) if (x+y) % 2 == 0)
    odd = ((1 << N) - 1) ^ even
    nine = sum(1 << (y*6+x) for y in (0, 2, 4) for x in (0, 2, 4))
    selected = {
        'corner': 1,
        'two_same_row': (1 << 0) | (1 << 2),
        'even_even_3x3': nine,
        'checkerboard_even': even,
        'checkerboard_odd': odd,
    }
    cpp = json.loads(subprocess.check_output([str(binary), '6'], text=True))
    got = {name: joint_count(mask) for name, mask in selected.items()}
    assert {name: str(value) for name, value in got.items()} == cpp['selected_joint_peak_counts']
    print(json.dumps({'status': 'passed', 'counts': {k: str(v) for k, v in got.items()}}, indent=2))


if __name__ == '__main__':
    main(Path(sys.argv[1]))
