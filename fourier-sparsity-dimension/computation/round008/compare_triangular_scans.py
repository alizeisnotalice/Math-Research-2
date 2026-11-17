#!/usr/bin/env python3
"""Compare all weight bins and all coordinate-activity bins across machines."""
import argparse
import hashlib
import json
from pathlib import Path

ap=argparse.ArgumentParser();ap.add_argument('local',type=Path);ap.add_argument('remote',type=Path);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args()
def aggregate(path):
    receipt=json.loads(path.read_text())
    assert receipt.get('execution_exit_code',receipt.get('ssh_exit_code'))==0
    counts=[0]*257;activity=[0]*257
    for shard in receipt['result']['shards']:
        assert shard['exit_code']==0
        d=shard['result']
        for w in range(257):counts[w]+=d['histogram'][w];activity[w]|=d['parameter_activity_by_weight'][w]
    assert sum(counts)==1<<31
    return counts,activity
local=aggregate(a.local);remote=aggregate(a.remote);assert local==remote
out={'status':'passed','codewords_per_machine':1<<31,'histogram_bins_compared':257,'activity_bins_compared':257,
     'local_receipt_sha256':hashlib.sha256(a.local.read_bytes()).hexdigest(),'remote_receipt_sha256':hashlib.sha256(a.remote.read_bytes()).hexdigest(),
     'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
a.out.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out))
