"""Compare C++ modular vectors to exact U192 and frozen compact-mask gates."""
import json
from pathlib import Path
import subprocess
import sys

from run_crt import HERE, MIB, crt, digest, independent_sizes, plan, verified_primes


def call(name, *args, timeout=45):
    proc = subprocess.run([str(HERE / name), *map(str, args)], check=True, capture_output=True,
                          text=True, timeout=timeout)
    assert proc.stderr == ''
    return json.loads(proc.stdout)


def exact_local_reference(n, h):
    out = call('u192_reference', n, h)
    assert out['profile'] == h and out['n'] == n
    assert out['independent_sets'] == sum(independent_sizes(h))
    return list(map(int, out['root_counts']))


def check_modular_vector(h, expect, primes):
    residues = []
    for ix, p in enumerate(primes):
        result = call('root_mod', h, ix)
        assert int(result['prime']) == p and result['prime_index'] == ix and result['profile'] == h
        histogram = independent_sizes(h)
        assert result['processed_by_layer'] == result['states_by_layer'] == histogram
        residue = list(map(int, result['root_residues']))
        assert residue == [a % p for a in expect]
        residues.append(residue)
    assert [crt(row, primes) for row in zip(*residues)] == expect
    return {'profile': h, 'n': len(h), 'independent_sets': sum(independent_sizes(h)),
            'root_counts': list(map(str, expect)), 'all_three_primes_and_reconstruction': True}


def main():
    primes, count = verified_primes(HERE / 'pratt-61bit.json')
    small = json.loads((HERE / 'reference-small.json').read_text())
    planned = json.loads((HERE / 'reference-plan.json').read_text())
    ranges = ['111', '122', '2323', '34343', '33434', '434343', '4545434', '44444444']
    gates = [dict(check_modular_vector(h, exact_local_reference(len(h), h), primes),
                  reference='fresh U192 independent-set implementation') for h in ranges]
    for row in small['small_readback']:
        if row['n'] > 8:
            gates.append(dict(check_modular_vector(row['profile'], list(map(int, row['rooted_counts'])), primes),
                              reference='frozen independent subset completion'))
    profiles_to_plan = {s['profile']: s for record in planned['profile_pairs'] for s in record['sides']}
    for h, row in profiles_to_plan.items():
        output = plan(h, 1024)
        assert row['independent_sets_by_size'] == output['states_by_layer']
        assert row['two_table_payload_bytes_64_bit_residue_plus_u64_key'] == output['max_table_bytes']
    larger = '345656543'
    constrained = plan(larger, 1023)
    assert constrained['max_table_bytes'] == 1024 * MIB and not constrained['fits_budget']
    reject = subprocess.run([str(HERE / 'root_mod'), larger, '0', '--max-table-mib', '1023'],
                            capture_output=True, text=True, timeout=10)
    assert reject.returncode == 2 and reject.stdout == '' and 'budget refused' in reject.stderr
    readback_path = HERE / 'n8-large-residue.json'
    reference_path = HERE / 'reference-n8-large.json'
    observed = json.loads(readback_path.read_text())
    reference = json.loads(reference_path.read_text())
    assert observed['profile'] == reference['profile']
    assert observed['independent_sets'] == reference['independent_sets'] == 10610055
    assert [int(a) % int(observed['prime']) for a in reference['root_counts']] == list(map(int, observed['root_residues']))
    small_batch = json.loads((HERE / 'small-driver-gate/summary.json').read_text())
    assert small_batch['completed_jobs'] == 9
    assert small_batch['root_counts'] == {g['profile']: g['root_counts'] for g in gates if g['profile'] in small_batch['root_counts']}
    audit = {'scope': 'C++ counting gates through n=8 and a sparse n=9 input, balanced n=9 planning only.',
             'python': sys.version, 'primes': list(map(str, primes)), 'pratt_nodes': count,
             'orig_small_gate_sha256': digest(HERE / 'reference-small.json'),
             'orig_n9_plan_sha256': digest(HERE / 'reference-plan.json'),
             'new_counter_sha256': digest(HERE / 'root_mod.cpp'),
             'binary_sha256': digest(HERE / 'root_mod'),
             'fresh_u192_reference_sha256': digest(HERE / 'u192_reference.cpp'),
             'verified_compact_n9_profile_plans': len(profiles_to_plan),
             'checked_1023_mib_refusal': True, 'verified_n8_largest_old_root_vector': True,
             'small_driver_completed_jobs': small_batch['completed_jobs'], 'counted_gates': gates}
    (HERE / 'gate-audit.json').write_text(json.dumps(audit, indent=2) + '\n')
    print(len(gates), 'gates with CRT; n9 plan', len(profiles_to_plan), 'regions; n8 large single residue verified')


if __name__ == '__main__':
    main()
