#!/usr/bin/env python3
"""Check valid ID 0 and three malformed replicate/RB ID pairs."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys

ap = argparse.ArgumentParser()
ap.add_argument('validator',type=Path)
ap.add_argument('raw',type=Path)
ap.add_argument('trace',type=Path)
ap.add_argument('source',type=Path)
ap.add_argument('output',type=Path)
args = ap.parse_args()
base = [json.loads(l) for l in args.raw.read_text().splitlines()]
records = []
for pair in ((0,0),(1,1),(0,1),(1,0)):
    rows = json.loads(json.dumps(base))
    rows[1]['replicate'],rows[2]['replicate'] = pair
    p = args.output.parent/f'probe-{pair[0]}-{pair[1]}.jsonl'
    assert not p.exists()
    p.write_text(''.join(json.dumps(r)+'\n' for r in rows))
    try:
        child = subprocess.Popen([sys.executable,'-B',str(args.validator),str(p),str(args.trace),str(args.source)],stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True)
        out,err = child.communicate()
        accepted = pair == (0,0)
        assert (child.returncode == 0) == accepted
        check = subprocess.run(['ps','-p',str(child.pid),'-o','pid='],capture_output=True,text=True)
        assert check.returncode == 1 and not check.stdout.strip()
        records.append({'ids':pair,'expected_accepted':accepted,'returncode':child.returncode,
            'stderr':err.strip(),'pid':child.pid,'exit_confirmed':True})
    finally:
        p.unlink()
result = {'status':'passed','scope':'trace ID validation only; production kernel unchanged',
    'validator_sha256':hashlib.sha256(args.validator.read_bytes()).hexdigest(),
    'test_source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    'raw_sha256':hashlib.sha256(args.raw.read_bytes()).hexdigest(),'tests':records}
args.output.write_text(json.dumps(result,indent=2)+'\n')
print('all four ID cases passed')
