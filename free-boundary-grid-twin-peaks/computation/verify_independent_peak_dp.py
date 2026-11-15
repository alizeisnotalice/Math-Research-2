"""Cross-check the independent-set peak DP against prior exact certificates.

The n=6 values are the newly produced certificate; n=2..5 are independent
checks against brute force or the pre-existing subset/peak-pair computations.
"""
import json
from pathlib import Path
import subprocess
import sys
from math import factorial
from fractions import Fraction


EXPECTED = {
    2: (16, 8, {2: 8}, 16),
    3: (29568, 143808, {2: 63552, 3: 64512, 4: 15744}, 633600),
    4: (35704394880, 779299174080,
        {2: 173341437440, 3: 313977112960, 4: 217620337600,
         5: 66746375360, 6: 7613910720}, 5272309842240),
    5: (88366736882654697600, 7063261481362958280000,
        {2: 759127145926367076000, 3: 2034099833145681772800,
         4: 2315729092924013995200, 5: 1393059798442968360000,
         6: 473582366731277968800, 7: 82379424872165332800,
         8: 5283819320483774400}, 67631454448823894520000),
    6: (1080582891996915763252912268789760,
        304351937651819767235871772704107520,
        {2: 14481687754950520311753606334274560,
         3: 54902115403376027519297614173153280,
         4: 89213621217858816803969061471462400,
         5: 80315330184021628648088627346048000,
         6: 45250034448130015901582999471764480,
         7: 16308868943054433429138573390192640,
         8: 3485524715479549508311417265587200,
         9: 379508572458565080597843375462400,
         10: 15246412490210033132029876162560},
        3934566741149040698531461804200320000),
}


def main(binary: Path):
    for n, (a1, a2, dist, sq) in EXPECTED.items():
        got = json.loads(subprocess.check_output([str(binary), str(n)], text=True))
        assert int(got['A1']) == a1
        assert int(got['A2']) == a2
        assert {int(k): int(v) for k, v in got['distance_counts'].items()} == dist
        assert int(got['squared_euclidean_total']) == sq
        histogram = {int(k): int(v) for k, v in got['peak_histogram'].items()}
        assert all(v > 0 for v in histogram.values())
        assert sum(histogram.values()) == factorial(n*n)
        assert sum(map(int, got['single_peak_root_counts'])) == a1
        assert sum(map(int, got['distance_counts'].values())) == a2
        ek = sum(Fraction(k*v, factorial(n*n)) for k, v in histogram.items())
        ek2 = sum(Fraction(k*(k-1)//2*v, factorial(n*n)) for k, v in histogram.items())
        assert ek == Fraction(3*n*n+3*n+2, 15)
        if n >= 4:
            assert ek2 == (Fraction(n**4, 50) + Fraction(n**3, 25)
                           - Fraction(11*n*n, 450) - Fraction(n, 15)
                           - Fraction(13, 315))
        if n == 6:
            archived = json.loads((Path(__file__).parent / 'results/independent-peaks-n6-r004.json').read_text())
            for key in ('A1', 'A2', 'peak_histogram', 'distance_counts',
                        'squared_euclidean_total', 'selected_joint_peak_counts', 'single_peak_root_counts'):
                assert archived[key] == got[key]
    print(json.dumps({'verified_n': list(EXPECTED), 'status': 'passed'}))


if __name__ == '__main__':
    main(Path(sys.argv[1]))
