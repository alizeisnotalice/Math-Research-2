"""Check archived/raw hashes and extract the explicitly listed ledgers to data temp."""
import argparse, gzip, hashlib, json, shutil
from pathlib import Path

HERE = Path(__file__).resolve().parent
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out',type=Path,required=True)
    args = ap.parse_args()
    if not str(args.out.resolve()).startswith('/Volumes/data/CodexTemp/'):
        raise SystemExit('Choose an output directory below mounted data/CodexTemp.')
    rows = json.loads((HERE/'compression-manifest.json').read_text())
    for row in rows:
        rel = Path(row['path'])
        assert rel.parent.name in {'search-d18-d19','search-d20-d24'}
        assert row['uncompressed_name'] in {'h-coverage.jsonl','signed-allocations.jsonl'}
        packed = (HERE/rel).read_bytes()
        assert hashlib.sha256(packed).hexdigest() == row['gzip_sha256']
        raw = gzip.decompress(packed)
        assert len(raw) == row['raw_bytes']
        assert hashlib.sha256(raw).hexdigest() == row['raw_sha256']
        dest = args.out/rel.parent
        dest.mkdir(parents=True,exist_ok=True)
        path = dest/row['uncompressed_name']
        if path.exists():
            assert path.read_bytes() == raw, 'Refusing to overwrite a different ledger'
        else:
            path.write_bytes(raw)
        for name in ['metadata.json','final.json']:
            data = (HERE/rel.parent/name).read_bytes()
            target = dest/name
            if target.exists():
                assert target.read_bytes() == data
            else:
                target.write_bytes(data)
    print({'verified_and_extracted_ledgers':len(rows),'out':str(args.out)})
if __name__ == '__main__':
    main()
