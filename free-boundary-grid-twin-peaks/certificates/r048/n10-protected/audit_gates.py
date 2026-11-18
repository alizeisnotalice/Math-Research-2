"""Gate a protected-peak four-prime build without a size-ten production run."""
import argparse
import json
from math import factorial, prod
from pathlib import Path
import subprocess

from count_protected import HERE, MIB, count, digest, frequencies, plan, verified_primes


def main():
    args = argparse.ArgumentParser()
    args.add_argument('old_ledger', type=Path)
    args.add_argument('old_pratt', type=Path)
    options = args.parse_args()
    primes, nodes = verified_primes(HERE/'pratt-61bit.json')
    assert primes[3] == 2305843009213693723 and nodes == 40
    assert prod(primes[:3]) < factorial(50) < prod(primes)
    verified = json.loads((HERE/'pratt-61bit.json').read_text())
    old_pratt = json.loads(options.old_pratt.read_text())
    assert verified['roots'][:3] == old_pratt['roots']
    assert all(verified['nodes'][q] == node for q,node in old_pratt['nodes'].items())
    ref = json.loads(options.old_ledger.read_text())
    replay_path = HERE/'small-replay'/'summary.json'
    replay = json.loads(replay_path.read_text())
    assert replay['completed_jobs'] == 16
    fields = ('counts', 'unsafe_max', 'unsafe_second_after_safe_nonadjacent_max')
    for h,vec in replay['region_histograms'].items():
        for source, expected in zip(vec, (ref['regions'][h][key] for key in fields)):
            expected += ['0'] * max(0, len(source)-len(expected))
            source += ['0'] * max(0, len(expected)-len(source))
            assert source == expected, h
    for row in replay['upper_counts']:
        old = next(z for z in ref['upper_bounds'] if z['profile'] == row['cut'])
        assert old['protected_top_two_refined_upper_A2'] == row['A2_upper_count']
        assert old['unprotected_nonadjacent_global_second_removed_words'] == row['global_second_excluded_count']
    straight = plan('5555555555', 512, primes)
    assert (straight['n'], straight['area'], straight['protected_vertices'], straight['independent_sets']) == (10, 50, 40, 32_641_916)
    assert straight['max_table_bytes'] == 512*MIB and max(straight['states_by_layer']) == 6_905_186
    refusal = subprocess.run([str(HERE/'protected_mod'), '5555555555', '3', '--max-table-mib', '511'],
                             capture_output=True, check=False, text=True, timeout=5)
    assert refusal.returncode == 2 and refusal.stdout == '' and 'two-table budget refused' in refusal.stderr
    # This profile has only one protected point, with two neighbors. The top
    # status counts follow from the rank symmetry of these three vertices.
    profile = '1111111112'
    sparse = HERE/'sparse-n10-gate'
    sparse.mkdir(exist_ok=False)
    p = plan(profile, 1, primes)
    parts = [count(profile,i,prime,sparse,p,1,10) for i,prime in enumerate(primes)]
    result = frequencies(profile, {k:r for k,r,_ in parts}, primes)
    m = sum(map(int, profile))
    yes_peak = factorial(m)//3
    yes_unsafe_top = (m-3)*factorial(m-1)//3
    yes_bad_second = (m-3)*factorial(m-2)
    expected = ([factorial(m)-yes_peak,yes_peak], [(m-1)*factorial(m-1)-yes_unsafe_top,yes_unsafe_top],
                [0,yes_bad_second])
    assert result == expected
    audit = {'scope':'n3/n9 four-modular replay plus sparse n10 test, n10 straight planning only',
             'prime_cert_sha256':digest(HERE/'pratt-61bit.json'),'verified_pratt_nodes':nodes,
             'previous_certificate_sha256':digest(options.old_pratt),'previous_ledger_sha256':digest(options.old_ledger),
             'prior_histogram_profiles_exactly_matched':sorted(replay['region_histograms']),
             'prior_upper_counts_exactly_matched':[f['cut'] for f in replay['upper_counts']],
             'actual_replay_sha256':digest(replay_path),'straight_plan':straight,
             'capacity_511_MiB_refusal_verified':True,'sparse_n10_profile':profile,
             'sparse_n10_histograms':result, 'sparse_n10_receipts':[rec for _,_,rec in parts]}
    (HERE/'gate-audit.json').write_text(json.dumps(audit,indent=2)+'\n')
    print('n3/n9 replay, 40 Pratt nodes, n10 plan and sparse gate passed')

if __name__ == '__main__':
    main()
