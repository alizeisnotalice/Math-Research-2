#!/usr/bin/env python3
"""Check the two archived n=6 thread schedules without remote inputs."""
import json
from pathlib import Path
here=Path(__file__).resolve().parent
a=json.loads((here/'n6-t1-check.json').read_text())
b=json.loads((here/'n6-t12.json').read_text())
fields=('n','M','algorithm','states_visited','H_two_rows','H_toprow','limit30')
assert {k:a[k] for k in fields}=={k:b[k] for k in fields}
assert a['n']==6 and a['M']==24 and a['states_visited']==16777216
assert a['H_two_rows']=='6637177983794054400'
assert a['H_toprow']=='20779014929894865426585600'
assert json.loads((here/'n6-t12.exit.json').read_text())['exit_code']==0
print(json.dumps({'passed':True,'comparison_fields':list(fields),'n6':{k:a[k] for k in fields}},indent=2))
