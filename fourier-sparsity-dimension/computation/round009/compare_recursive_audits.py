#!/usr/bin/env python3
"""Compare all exact certificates across local and remote executions."""
import argparse,hashlib,json
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);a=p.parse_args();base=Path(__file__).parent
local_path=base/'local-recursive-audit.json';remote_path=base/'remote-recursive-audit.json'
local=json.loads(local_path.read_text());receipt=json.loads(remote_path.read_text())
assert receipt['ssh_exit_code']==receipt['remote']['exit_code']==0
remote=receipt['remote']['result'];assert local['status']==remote['status']=='passed'
source_hash=hashlib.sha256((base/'recursive_normal_rank_certificate.py').read_bytes()).hexdigest()
assert local['source_sha256']==remote['source_sha256']==receipt['remote']['source_sha256']==source_hash
sections=['n4_c_le_2_exhaustive','mixed_depth_tree_trials','full_tree_equality_controls']
for key in sections:assert local[key]==remote[key],key
left=dict(local['r008_remote_packing_witnesses']);right=dict(remote['r008_remote_packing_witnesses'])
left.pop('receipt_path');right.pop('receipt_path');assert left==right
assert left['receipt_sha256']==receipt['remote']['input_sha256']==hashlib.sha256((base.parent/'round008/remote-packing-execution.json').read_bytes()).hexdigest()
result={'status':'passed','comparison':'all exact certificate fields; only host Python/time and input pathname excluded',
 'exhaustive_families':local['n4_c_le_2_exhaustive']['families_including_empty'],
 'exhaustive_certificate_digest':local['n4_c_le_2_exhaustive']['case_certificate_rolling_sha256'],
 'tree_subfamilies':local['mixed_depth_tree_trials']['trial_count'],'equality_controls':len(local['full_tree_equality_controls']['controls']),
 'previous_packing_witnesses':len(left['records']),'local_seconds':local['elapsed_seconds'],
 'remote_seconds':remote['elapsed_seconds'],'remote_host':receipt['remote']['host'],
 'local_input_sha256':hashlib.sha256(local_path.read_bytes()).hexdigest(),
 'remote_input_sha256':hashlib.sha256(remote_path.read_bytes()).hexdigest(),
 'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
a.out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
