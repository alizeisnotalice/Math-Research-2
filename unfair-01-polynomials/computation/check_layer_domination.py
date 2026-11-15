#!/usr/bin/env python3
"""Check frozen domination certificates and exact coverage of all unknowns."""
import argparse
import json
from pathlib import Path
from diagonal_domination import verify_certificate


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('ledger_metadata',type=Path)
    ap.add_argument('certificates',type=Path)
    args=ap.parse_args()
    metadata=json.loads(args.ledger_metadata.read_text())
    cases=json.loads(args.certificates.read_text())['cases']
    key=lambda c:(tuple(c['I']),tuple(c['J']))
    unknown=[key(c) for g in metadata['groups'] for j in g['jobs'] for c in j['unknown']]
    assert len(set(unknown))==len(unknown)
    assert sorted(map(key,cases))==sorted(unknown)
    for c in cases:
        assert c['status']=='excluded'
        verify_certificate(c['I'],c['J'],c['certificate'])
    print(json.dumps({'all_unknowns_exactly_covered':True,'verified_certificates':len(cases)}))


if __name__=='__main__':main()
