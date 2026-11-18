"""Run a single independently planned modulus and save lifecycle receipts."""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path

from count_protected import count, digest, plan, verified_primes

def utc():
    return datetime.now(timezone.utc).isoformat()

def main():
    p = argparse.ArgumentParser()
    p.add_argument('profile')
    p.add_argument('output', type=Path)
    p.add_argument('--max-table-mib', required=True, type=int)
    p.add_argument('--timeout', type=int, default=10800)
    p.add_argument('--prime-index', type=int, default=0)
    args = p.parse_args()
    primes, nodes = verified_primes(Path('pratt-61bit.json'))
    assert nodes == 50
    idx = args.prime_index
    assert 0 <= idx < len(primes)
    checked = plan(args.profile,args.max_table_mib,primes)
    assert checked['fits_budget'] and checked['five_moduli_cover_factorial']
    args.output.mkdir()
    meta = dict(started_utc=utc(), source_sha256=digest(Path('protected_mod.cpp')),
                binary_sha256=digest(Path('protected_mod')), prime_nodes=nodes, prime_index=idx,
                timeout=args.timeout, planned=checked)
    (args.output/'started.json').write_text(json.dumps(meta,indent=2)+'\n')
    print('accepted',args.profile,idx,checked['max_table_bytes'],flush=True)
    try:
        _, _, receipt = count(args.profile,idx,primes[idx],args.output,checked,
                              args.max_table_mib,args.timeout)
        meta.update(receipt=receipt,outcome='complete')
    except BaseException as e:
        meta.update(outcome='error',error=str(e))
        raise
    finally:
        meta.update(ended_utc=utc())
        (args.output/'finished.json').write_text(json.dumps(meta,indent=2)+'\n')
    print('finished and waited',args.profile,idx,receipt['max_rss_bytes'],flush=True)

if __name__ == '__main__':
    main()
