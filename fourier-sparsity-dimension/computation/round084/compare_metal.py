#!/usr/bin/env python3
"""Compare the signed integer GPU data, ranks, and sparsities across both hosts."""
import json
from pathlib import Path


def load(name):
    return json.loads(Path(name).read_text())


def main():
    out=[]
    fields=('seed','mode','stage','initial_k','coefficient_count','s','D','P_squared','actual_int32_walsh_sha256')
    for prefix in ('phase-metal','extras-metal'):
        local,remote=[load(f'{prefix}-{host}.json') for host in ('mini','macbook')]
        left,right=[[[r[f] for f in fields] for r in data['records']] for data in (local,remote)]
        if left!=right:
            raise AssertionError(('different GPU SHA-256 signatures or counts',prefix))
        out.append(dict(batch=prefix,exact_lifts=len(left),coefficients=sum(x[4] for x in left),
                        devices=[z['device'] for z in (local,remote)],
                        gpu_command_seconds=[z['gpu_command_seconds'] for z in (local,remote)],
                        signed_spectrum_sha256=[r[-1] for r in left]))
    print(json.dumps(dict(status='pass',identical_batches=out,
                          coefficients_per_device=sum(x['coefficients'] for x in out)),indent=2,sort_keys=True))


if __name__=='__main__':
    main()
