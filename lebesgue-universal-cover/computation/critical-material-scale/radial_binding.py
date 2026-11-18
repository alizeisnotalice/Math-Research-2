from pathlib import Path
from hashlib import sha256
import numpy as np,json,sys
c=Path(sys.argv[1]).resolve();t=Path(sys.argv[2]).resolve();out=Path(sys.argv[3]);rec=json.loads((c/'five-body-radial-affine.json').read_text());rad=c/'five-body-radial-material.npz';th=sha256(t.read_bytes()).hexdigest();rh=sha256(rad.read_bytes()).hexdigest();assert th==rec['source_table_sha256'] and rh==rec['table_sha']
with np.load(t,allow_pickle=False) as a,np.load(rad,allow_pickle=False) as b:
 assert np.array_equal(a['counts'],b['counts']) and np.array_equal(a['points'][0],b['points'])
out.write_text(json.dumps({'source_table_sha256':th,'radial_table_sha256':rh,'zero_angle_and_counts_equal':True},indent=2)+'\n')
