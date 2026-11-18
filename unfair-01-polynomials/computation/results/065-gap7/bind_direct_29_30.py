"""Bind the two separately generated exact prefix files to the direct modular census."""
from argparse import ArgumentParser
from hashlib import file_digest
import json
from pathlib import Path


def sha(path: Path) -> str:
    with path.open('rb') as source:
        return file_digest(source, 'sha256').hexdigest()


def main():
    if not __debug__:
        raise SystemExit('requires assertions')
    parser=ArgumentParser(description=__doc__)
    parser.add_argument('--prefix', required=True, type=Path)
    parser.add_argument('--fresh-prefix', required=True, type=Path)
    parser.add_argument('--direct', required=True, type=Path)
    parser.add_argument('--arithmetic-audit', required=True, type=Path)
    parser.add_argument('--out', required=True, type=Path)
    args=parser.parse_args()
    prefix_hash=sha(args.prefix)
    assert prefix_hash==sha(args.fresh_prefix)
    direct=json.loads(args.direct.read_text())
    arithmetic=json.loads(args.arithmetic_audit.read_text())
    prefix=json.loads(args.prefix.read_text())
    assert arithmetic['status']=='verified' and arithmetic['input_sha256']==prefix_hash
    assert direct['mode']=='scout' and direct['status']=='complete' and direct['offset']==0
    assert [row['degree'] for row in direct['degrees']]==[29,30]
    for cases, result in zip(prefix,direct['degrees'],strict=True):
        assert result['degree']==cases['degree'] and cases['counts']==result['early_counts']
        assert result['modular_counts']['shapes']==len(cases['late'])
        assert not result['unresolved']
    root=Path(__file__).resolve().parent
    for name,digest in direct['source_sha256'].items():
        assert sha(root/name)==digest
    status={'status':'verified','original_and_fresh_prefix_sha256':prefix_hash,
            'direct_sha256':sha(args.direct),'arithmetic_audit_sha256':sha(args.arithmetic_audit),
            'source_sha256':direct['source_sha256'],'degrees':direct['degrees'],
            'direct_seconds':direct['seconds'],'python':direct['python'],'sympy':direct['sympy']}
    with args.out.open('x') as out:
        out.write(json.dumps(status,indent=2,sort_keys=True)+'\n')


if __name__=='__main__':
    main()
