#!/usr/bin/env python3
"""Exact K=2 peak-indicator pair-association audit on small square grids."""
import itertools
import json
import subprocess
import sys


def d1(n, p, q):
    return abs(p % n - q % n) + abs(p // n - q // n)


EXPECTED = {3: 143808, 4: 779299174080}


def run(binary, n):
    K = {}
    for p, q in itertools.combinations(range(n*n), 2):
        if d1(n, p, q) > 1:
            K[p, q] = int(subprocess.check_output(
                [binary, str(n), str(p), str(q)], text=True))
    A = sum(K.values()); assert A == EXPECTED[n]
    marg = [0]*(n*n)
    for (p, q), z in K.items():
        marg[p] += z; marg[q] += z
    # Pairwise negative association of peak indicators would say
    # K(p,q)/A <= (m_p/A)(m_q/A), i.e. residual <= 0.
    max_pos = None; min_res = None; positives = 0
    for (p, q), z in K.items():
        r = z*A - marg[p]*marg[q]
        item = {'p': [p % n, p // n], 'q': [q % n, q // n],
                'distance': d1(n, p, q), 'K': z,
                'm_p': marg[p], 'm_q': marg[q], 'residual': r}
        if r > 0:
            positives += 1
            if max_pos is None or r > max_pos['residual']:
                max_pos = item
        if min_res is None or r < min_res['residual']:
            min_res = item
    # Layer-cake cut signs: 2x-(n-1)=sum_{j=0}^{n-2}s_j(x), where s_j is
    # +1 right of cut j and -1 on/left of it.  Means need not be zero for
    # off-centre cuts, so audit the *centered* covariance matrix explicitly.
    cuts = []
    signs = [[1 if v % n > j else -1 for v in range(n*n)]
             for j in range(n-1)]
    means = [sum(marg[v]*s[v] for v in range(n*n)) for s in signs]
    for j in range(n-1):
        raw = sum(z*signs[j][p]*signs[j][q] for (p, q), z in K.items())
        cuts.append({'cut_after_x': j, 'mean_numerator_over_2A': means[j],
                     'raw_same_cut_numerator_over_A': raw})
    covariance = []
    for j in range(n-1):
        row = []
        for k in range(n-1):
            ordered_num = sum(z*(signs[j][p]*signs[k][q] +
                                 signs[j][q]*signs[k][p])
                              for (p, q), z in K.items())
            # Cov = ordered_num/(2A)-means[j]means[k]/(4A^2).
            centered_num = 2*A*ordered_num - means[j]*means[k]
            row.append({'numerator_over_4A2': centered_num,
                        'sign': ('positive' if centered_num > 0 else
                                 ('negative' if centered_num < 0 else 'zero'))})
        covariance.append(row)
    return {'n': n, 'A2': A, 'nonadjacent_pairs': len(K),
            'positive_pair_association_count': positives,
            'largest_positive_association': max_pos,
            'most_negative_association': min_res,
            'cut_moments': cuts,
            'centered_cut_covariances': covariance,
            'marginals': marg}


if __name__ == '__main__':
    if len(sys.argv) != 3:
        raise SystemExit('usage: audit_peak_association.py PAIR_DP_BINARY OUT')
    out = {'results': [run(sys.argv[1], n) for n in (3, 4)]}
    with open(sys.argv[2], 'w') as f:
        json.dump(out, f, indent=2); f.write('\n')
    print(json.dumps(out, indent=2))
