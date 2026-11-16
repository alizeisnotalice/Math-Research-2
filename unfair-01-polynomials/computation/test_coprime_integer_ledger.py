"""Reject five malformed archived-ledger fixtures without invoking the GPU."""
import argparse
import copy
import json
from pathlib import Path
from audit_coprime_integer_ledger import audit, sha

ap=argparse.ArgumentParser()
ap.add_argument('--archive-dir',type=Path,required=True)
args=ap.parse_args()
root=args.archive_dir
c=json.loads((root/'checkpoints/a2b9-checkpoint.json').read_text())
raw=json.loads((root/'inputs/a2b9-input.json').read_text())['candidates'][0]
source=sha(root/'scripts/stream_coprime_ds_gpu_v2.py')
mutations={
    'duplicate_saved_code':lambda x:x['E_progress'][0]['top64'].__setitem__(1,copy.deepcopy(x['E_progress'][0]['top64'][0])),
    'wrong_count':lambda x:x['E_progress'][0]['counts'].__setitem__('codes_seen',0),
    'wrong_mask':lambda x:x['E_progress'][0]['top64'][0].__setitem__('A_root_group_bitmask_decimal','0'),
    'missing_E':lambda x:x['E_progress'].pop(0),
}
for name,mutate in mutations.items():
    changed=copy.deepcopy(c)
    mutate(changed)
    try: audit(changed,raw,source)
    except (AssertionError,ValueError,KeyError): pass
    else: raise AssertionError(f'Invalid fixture accepted: {name}')
changed=copy.deepcopy(raw)
changed['blocks'][1]['a_lift_groups'][0]=changed['blocks'][0]['a_lift_groups'][0]
try: audit(c,changed,source)
except (AssertionError,ValueError,KeyError): pass
else: raise AssertionError('Duplicate root index accepted')
print(json.dumps({'status':'all_five_rejected','cases':list(mutations)+['duplicate_root_index']}))
