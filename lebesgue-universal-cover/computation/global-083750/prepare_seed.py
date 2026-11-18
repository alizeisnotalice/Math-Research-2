import argparse,hashlib,json,sys
from pathlib import Path
from fractions import Fraction as F
B,C=map(Path,sys.argv[1:]);sys.path.insert(0,str(C.resolve()))
from prepare_tomar_multi_target_layer import prepare
base=json.loads((B/'base-report.json').read_bytes());mr=json.loads((B/'merged45-report.json').read_bytes());assert mr['remaining']==699060
can={r[0]:r[1] for r in base['rows']};base['rows']=[[p,can[p],hp,hf] for p,n,hp,hf in mr['rows']];base['curve']=[];base['provenance']['frozen_score_merge']=hashlib.sha256((B/'merged45-report.json').read_bytes()).hexdigest();report=B/'full-index-report.json';report.write_text(json.dumps(base,separators=(',',':'))+'\n')
r=prepare(argparse.Namespace(index_report=report,index=B/'merged45-index',radial_certificate=C/'five-body-radial-affine.json',emit=B/'seeds',target=F('.83750'),radial_target=F('.83750')));assert r['remaining']==699060;(B/'seed.json').write_text(json.dumps(r,separators=(',',':'))+'\n');print('remaining',r['remaining'])
