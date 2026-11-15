#!/usr/bin/env python3
"""Hash current R61 artifacts for root handoff; independent review remains pending.

Prepared but not run by the author after the permission-profile change left its
executor unable to start. Root may run this in its refreshed full-access session.
"""
import hashlib
import json
import os
from pathlib import Path

assert os.path.ismount('/Volumes/data')
ROOT = Path('/Volumes/data/CodexProjects/T2')
WORKSPACES = ROOT / '.local/workspaces'
TOPIC = 'planar-centered-square-weak11'
OWN = WORKSPACES / 'pcw-astra-recursive' / TOPIC
OPS = WORKSPACES / 'pcw-operations-terra/r061'
own_paths = [
    'proofs/finite-source-union-cut-certificate.md',
    'rounds/061.md',
    'computation/round061/finite_source_cut_certificate.py',
    'computation/round061/r061-bounded-inputs.json',
    'computation/round061/summarize_finite_source_cuts.py',
    'computation/round061/finite-source-cut-summary.json',
    'computation/round061/critical-input-minimal.json',
    'computation/round061/finalize_round061_manifest.py',
    'reviews/round061-future-frontier-handoff.md',
]
files = []


def add(source, destination, expected=None):
    data = source.read_bytes()
    digest = hashlib.sha256(data).hexdigest()
    assert len(data) < 10 * 1024 * 1024, source
    if expected is not None:
        assert digest == expected, (source, digest, expected)
    files.append({'source': str(source), 'destination': destination,
                  'sha256': digest, 'bytes': len(data)})


pins = {
    'computation/round061/finite_source_cut_certificate.py':
        'e1d356ac0771b36e3a2145775e98ee27d3211c56447f29fca6939e2b540e9b41',
    'computation/round061/r061-bounded-inputs.json':
        'c583b0a80014dca0b4e0d53b74f26bde2a98c9f335e9adbdc14315473749def2',
    'computation/round061/summarize_finite_source_cuts.py':
        'b2e49598d65b6480385b03b3d13016c660017033588cb09d54dc7a8477493d50',
}
for relative in own_paths:
    add(OWN / relative, TOPIC + '/' + relative, pins.get(relative))

summary = json.loads((OWN / 'computation/round061/finite-source-cut-summary.json').read_text())
source_pins = {Path(row['path']).name: row['sha256'] for row in summary['batch_sources']}
for batch in range(12):
    name = f'batch-{batch}.json'
    add(OPS / name, TOPIC + '/computation/results/round061/' + name,
        source_pins[name])
add(OPS / 'run-receipt.json', TOPIC + '/computation/results/round061/run-receipt.json',
    '7b1b0d100cec7364aabc22b23ed0f82f4b4058b52d412d011d024cdf2ba78525')
add(OPS / 'run_batches.py', TOPIC + '/computation/round061/run_batches.py')

sources = []
for workspace, relative, digest in [
    ('pcw-sol-cover-counterexample',
     'proofs/r053-nine-cell-mass-simplex-certificate-33over10.md',
     'f7e9223403939788d940ebc1cfd47d69ad772eb44da73e42cdeb08d6f80fd784'),
    ('pcw-astra-covering',
     'reviews/r053-covering-core-independent-review.md',
     'f466f85312c98c116f8276a10b14fbde3f8d9870be1ac5f970c4e196d76ece66'),
]:
    source = WORKSPACES / workspace / TOPIC / relative
    actual = hashlib.sha256(source.read_bytes()).hexdigest()
    assert actual == digest, (source, actual)
    sources.append({'source': str(source), 'existing_destination': TOPIC + '/' + relative,
                    'sha256': digest, 'scope': 'R53 geometric theorem only; not rerun in R61'})

assert len({f['destination'] for f in files}) == len(files)
manifest = {
    'round': 61,
    'date': '2026-11-15',
    'status': 'AUTHOR-CHECKED; INDEPENDENT REVIEW PENDING',
    'scope': 'Positive finite-input source-union cut and legal nonproportional allocation certificates; no new general C2 bound',
    'files': files,
    'total_payload_bytes': sum(f['bytes'] for f in files),
    'existing_mathematical_sources': sources,
    'actual_runs': {
        'primary': '12 Python 3.13.18 processes, all exit 0; coordinator confirmed all owned PIDs exited',
        'author_summary': 'Python 3.15.0, exact Fraction/integer Gamma and duplicate accounting, exit 0',
        'cases_evaluated': 768,
        'distinct_inputs': 725,
        'repeated_controls': 43,
        'nonempty_cuts': 1217024,
        'compact_cuts': 397434,
        'noncompact_cuts': 819590,
        'noncompact_critical_cases': 165,
        'terminal_improved_cuts': 0,
    },
    'independent_review': {
        'completed': False,
        'reason': 'Messages to completed R56 helpers did not start R61 work. One bounded restart was interrupted immediately under the new user team cap. No completed binding review or independent audit result is claimed.',
        'next_action': 'Selected next-phase team should review proof, geometry sweep, and all saved cut/flow certificates before publication.',
    },
    'future_frontier': '16–20-source plan saved, not implemented and not run',
    'git': 'No delegated commit or push; root retains integration, review and publication responsibility',
}
output = OWN / 'computation/results/round061-integration-manifest.json'
output.parent.mkdir(parents=True, exist_ok=True)
output.write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + '\n')
print(json.dumps({'path': str(output), 'files': len(files),
                  'payload_bytes': manifest['total_payload_bytes'],
                  'sha256': hashlib.sha256(output.read_bytes()).hexdigest(),
                  'review_status': manifest['status']}, indent=2))
