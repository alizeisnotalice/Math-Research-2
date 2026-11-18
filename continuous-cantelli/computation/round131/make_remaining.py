"""Combine adjacent unresolved rectangle cells for the exact Simpson stage."""
from decimal import Decimal
import json
from pathlib import Path

out={}
for path in Path('broad').glob('shard-*.json'):
    for row in json.loads(path.read_text())['results']:
        spans=[]
        for a,b,_ in sorted(row['unresolved'],key=lambda t:Decimal(t[0])):
            if spans and Decimal(spans[-1][-1])==Decimal(a):
                spans[-1][-1]=b
            else:
                spans.append([a,b])
        if spans:
            out[row['j']]=spans
Path('remaining.json').write_text(json.dumps(out,indent=2)+'\n')
