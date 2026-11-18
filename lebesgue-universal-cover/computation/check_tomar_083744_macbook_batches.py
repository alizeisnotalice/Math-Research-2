"""Check the six MacBook arithmetic runs against the complete 0.83744 layer ledger."""
from hashlib import sha256
import json
from pathlib import Path
from zipfile import ZipFile

HERE=Path(__file__).resolve().parent


def main():
    if not __debug__:
        raise SystemExit('Assertions must be enabled.')
    sha=lambda b:sha256(b).hexdigest()
    ledger=json.loads((HERE/'tomar-083744-stage-ledger.json').read_bytes())
    bindings=(HERE/'tomar-083744-macbook-batch-bindings.json').read_bytes()
    assert ledger['binary_sha']==sha(bindings)
    parts=set()
    aggregated=[{} for _ in range(4)]
    members=[[] for _ in range(4)]
    pins=json.loads(bindings)
    with ZipFile(HERE/'tomar-083744-macbook-batches.zip') as z:
        assert len(z.namelist())==18 and len(pins)==6
        for key,b in sorted(pins.items()):
            i=key.rsplit('-',1)[1]
            docs={}
            for field,name in [('geometry_ledger','stage-ledger.json'),('result_sha256','result.json'),('group_sha256','group.json')]:
                raw=z.read(f'batch-{i}/'+name)
                assert sha(raw)==b[field]
                docs[name]=json.loads(raw)
            g=docs['group.json']['parts']
            assert len(g)==len(set(g))==b['part_count'] and parts.isdisjoint(g)
            parts.update(g)
            stage=docs['stage-ledger.json']
            assert stage['binary_sha']==b['binary_sha256'] and stage['parts']==len(g)
            assert all(stage[k]==ledger[k] for k in ('radial_target','integer_goal','seed_receipt_sha','source_sha','target'))
            result=docs['result.json']
            assert result['complete'] and result['ledger_sha256']==b['geometry_ledger']
            assert result['group_sha256']==b['group_sha256'] and result['parts']==len(g)
            for j,a in enumerate(stage['audits']):
                assert [p for p,h in a['work_files']]==g
                members[j].extend(a['work_files'])
                for k,v in a.items():
                    if k!='work_files':
                        aggregated[j][k]=aggregated[j].get(k,0)+v
    assert len(parts)==ledger['parts']==1024
    for j,stage in enumerate(ledger['audits']):
        assert sorted(members[j])==stage['work_files']
        assert aggregated[j]=={k:v for k,v in stage.items() if k!='work_files'}
    print(json.dumps(dict(complete=True,arithmetic_runs=len(pins),parts=len(parts),
                          calls=sum(a['calls'] for a in aggregated),remaining=aggregated[-1]['remaining'])))


if __name__=='__main__':
    main()
