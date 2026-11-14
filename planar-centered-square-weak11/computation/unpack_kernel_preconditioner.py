"""Restore the exact 1182-square preconditioner, without changing the checker.

python unpack_kernel_preconditioner.py ARCHIVE.byteshuffle8.gz OUTPUT.bin
The storage manifest must be beside ARCHIVE. OUTPUT should be in CodexTemp.
"""
import gzip
import hashlib
import json
from pathlib import Path
import sys

archive,output=map(Path,sys.argv[1:])
meta=json.loads((archive.parent/'floor-preconditioner-storage.json').read_text())
packed=archive.read_bytes()
assert hashlib.sha256(packed).hexdigest()==meta['gzip_sha256']
planes=gzip.decompress(packed)
assert len(planes)==meta['raw_bytes'] and len(planes)%8==0
n=len(planes)//8
raw=bytearray(len(planes))
for k in range(8):
    raw[k::8]=planes[k*n:(k+1)*n]
assert hashlib.sha256(raw).hexdigest()==meta['raw_sha256']
output.parent.mkdir(parents=True,exist_ok=True)
output.write_bytes(raw)
print(json.dumps({'verified':True,'bytes':len(raw),'sha256':meta['raw_sha256'],'output':str(output)},indent=2))
