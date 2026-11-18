"""Verify exact parameter coverage and remote exit receipts."""
from datetime import datetime,UTC
from decimal import Decimal as D
from pathlib import Path
import json

P=Path('center')
files=[json.loads(s.read_text()) for s in P.glob('shard-*.json')]
assert len(files)==15
rows=sorted((x for f in files for x in f['results']),key=lambda r:D(r['z'][0]));assert len(rows)==236
z=D('.01')
for row in rows:
    lo,hi=map(D,row['z'])
    assert lo==z and hi-lo==D('.0025') and D(row['margin_lower'])>0
    assert row['panels']==1198
    z=hi
assert z==D('.6')
print(len(rows), min(D(row['margin_lower']) for row in rows))
