"""Combine the exact protected-count bound with the disjoint n9 numerator."""
from argparse import ArgumentParser
from fractions import Fraction
import json
from math import factorial
from pathlib import Path

from count_upper import decimal, digest


def main():
    p = ArgumentParser()
    p.add_argument('denominator', type=Path)
    p.add_argument('numerator', type=Path)
    p.add_argument('output', type=Path)
    arguments = p.parse_args()
    audit = json.loads(arguments.denominator.read_text())
    family = json.loads(arguments.numerator.read_text())
    candidates = [record for record in audit['upper_bounds'] if record['n'] == 9]
    chosen = min(candidates, key=lambda record: int(record['protected_top_two_refined_upper_A2']))
    upper = int(chosen['protected_top_two_refined_upper_A2'])
    strata = family['rootwise_counts']
    assert len({tuple(row['peaks']) for row in strata}) == len(strata) == family['root_pairs']
    total_by_dist = {}
    for row in strata:
        u, v = row['peaks']
        assert 0 <= u < v < 81
        distance = abs(u%9-v%9) + abs(u//9-v//9)
        assert distance > 9/3
        total_by_dist[distance] = total_by_dist.get(distance, 0) + row['count']
    eligible = [(Fraction(9, 3), sum(total_by_dist.values())),
                (Fraction(9, 2), sum(count for d,count in total_by_dist.items() if 2*d > 9))]
    assert eligible[0][1] == family['rootwise_subunion_count'] <= upper
    out = {'n': 9, 'profiles': len(family['profiles']), 'denominator_bound_selected_profile': chosen['profile'],
           'strict_A2_upper_count': str(upper), 'denominator_upper_over_81_factorial': decimal(Fraction(upper, factorial(81))),
           'denominator_provenance_sha256': digest(arguments.denominator),
           'numerator_provenance_sha256': digest(arguments.numerator),
           'rootwise_distance_counts': {str(d):str(count) for d,count in sorted(total_by_dist.items())},
           'conditional_lower_bound': [{'distance_strictly_larger_than': str(threshold), 'counted_words': str(count),
                                        'exact': str(Fraction(count, upper)), 'decimal': decimal(Fraction(count, upper))}
                                       for threshold, count in eligible]}
    with arguments.output.open('x') as dest:
        json.dump(out, dest, indent=2)
        dest.write('\n')
    for row in out['conditional_lower_bound']:
        print(row['distance_strictly_larger_than'], row['decimal'])


if __name__ == '__main__':
    main()
