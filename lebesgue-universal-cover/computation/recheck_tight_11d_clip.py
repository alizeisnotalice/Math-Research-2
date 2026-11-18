"""Check every metadata root and run the dyadic 11D clipping counter on frozen files."""
from argparse import ArgumentParser
from fractions import Fraction as Q
from hashlib import sha256
import json
from pathlib import Path
import subprocess

if __name__=='__main__':
    ap=ArgumentParser();ap.add_argument('directory',type=Path);ap.add_argument('executable',type=Path);a=ap.parse_args()
    source=Path(__file__).with_name('count_tight_11d_clip.cpp');baseline=None;hashes=sha256()
    for index in range(2048):
        m=json.loads((a.directory/f'part{index:04d}.json').read_bytes());assert m['part']==index and m['nparts']==2048 and m['split0']==[2]*11 and m['ns']==[3,5,7,9]
        root=(m['lo'],m['hi']);baseline=baseline or root;assert baseline==root
        hashes.update((a.directory/f'part{index:04d}.npy').read_bytes())
    numbers=[Q(x)*2**54 for seq in baseline for x in seq];assert all(x.denominator==1 for x in numbers)
    c=json.loads(subprocess.check_output([a.executable,a.directory,*map(str,numbers)]))
    assert c['leaves']==16193402
    print(json.dumps({'counts':c,'tree_stream_sha256':hashes.hexdigest(),'source_sha256':sha256(source.read_bytes()).hexdigest(),
                'angle_upper':'22/105','radii':['.0700','.0670','.0606','.0575'],'scope':'positive retention after conservative sphere and angle tests; does not verify area of retained leaves'},indent=2))
