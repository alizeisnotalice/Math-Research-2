#!/usr/bin/env python3
"""Focused semantic regressions for stream validation, no GPU run required."""
import copy, importlib.util, json
from pathlib import Path

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('stream',HERE/'native_stream_scheduler.py');mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
doc=json.loads(Path('/Volumes/data/CodexTemp/T2/unfair-01-polynomials/round038-macbook-gpu52-direct/inputs.json').read_text());task=doc['tasks'][2]
dims={'groups':task['groups'],'nfreq':len(doc['DC']),'ncoef':len(doc['DC'][0])};sha='a'*64;expected=[(0,16)]
base={'event':'batch','command_status':'completed','task':'2','input_sha256':sha,'source_sha256':sha,'binary_sha256':sha,'groups':dims['groups'],'nfreq':dims['nfreq'],'ncoef':dims['ncoef'],'range_start':0,'range_count':16,'finite_count':16,'top_indices':[0,1,2,3],'top_scores':[1.,1.,1.,1.]}
cases={
 'three_indices':lambda x:x.update(top_indices=[0,1,2]),
 'five_indices':lambda x:x.update(top_indices=[0,1,2,3,4]),
 'duplicate_indices':lambda x:x.update(top_indices=[0,1,1,3]),
 'fractional_index':lambda x:x.update(top_indices=[0,1.0,2,3]),
 'out_of_range_index':lambda x:x.update(top_indices=[0,1,2,16]),
}
passed=[]
for name,edit in cases.items():
 x=copy.deepcopy(base);edit(x)
 try:mod.validate_events([json.dumps(x)],expected,task,dims,sha,sha,sha)
 except ValueError:passed.append(name)
 else:raise AssertionError(name+' unexpectedly accepted')
bad_target=copy.deepcopy(task);bad_target['C'][0]=0
try:mod.direct_score(bad_target,1,dims['ncoef'])
except ValueError:passed.append('bad_target_dimension_or_nonfinite_rejected')
else:
 # a target mismatch must remain an explicit unknown, not acceptance
 assert mod.direct_score(bad_target,1,dims['ncoef'])['target_product_status']=='unknown_target_mismatch';passed.append('target_mismatch_unknown')
out={'status':'passed','negative_cases':passed,'count':len(passed)}
(HERE/'semantic_negative_tests.json').write_text(json.dumps(out,sort_keys=True,indent=2)+'\n');print(json.dumps(out,sort_keys=True))
