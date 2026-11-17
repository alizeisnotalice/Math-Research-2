from pathlib import Path
import argparse,gzip,hashlib,json,os
ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True);a=ap.parse_args()
assert os.path.ismount('/Volumes/data')
r=Path(__file__).resolve().parent;manifest=json.loads((r/'stdout-manifest.json').read_text());digest=hashlib.sha256();lines=0
with a.out.open('wb') as out:
 for part in manifest['parts']:
  packed=(r/part['path']).read_bytes();assert hashlib.sha256(packed).hexdigest()==part['gzip_sha256']
  raw=gzip.decompress(packed);assert hashlib.sha256(raw).hexdigest()==part['uncompressed_sha256'];assert lines==part['first_line']
  assert len(raw.splitlines())==part['line_count'];lines+=part['line_count'];digest.update(raw);out.write(raw)
assert lines==manifest['line_count'] and digest.hexdigest()==manifest['original_sha256']
print('verified lossless saved-stream reconstruction',lines,digest.hexdigest())
