"""Replay a constant-width pentagon and prove monotonic exclusion of larger pivots."""
from pathlib import Path
from fractions import Fraction as F
from hashlib import sha256
from flint import arb
import argparse,json,sys
from certify_sigma_perturbed_pentagon_obstruction import run
from certify_sigma_extended_graph import geometry

def main():
 ap=argparse.ArgumentParser();ap.add_argument('archive',type=Path);ap.add_argument('output',type=Path);a=ap.parse_args();comp=a.archive.resolve();cfgpath=comp/'pentagon-deep-pivot-parameters.json';cfg=json.loads(cfgpath.read_bytes());assert cfg=={'sigma':'1.5497','w':'.77347','deformation':'.03012','pivot':'.83331','angle':'.01262','cut_only':True}
 raw=a.output.with_suffix('.body.json');run(comp,cfgpath,raw);doc=json.loads(raw.read_bytes());ns,_,_=geometry(comp,[(F('.429'),F('.009'))],sigma=cfg['sigma'],far=F('.273'),negative=F('.306'),eps=F('.0102'));A=ns['A'];angle=A(F(cfg['angle']));direction=ns['add'](ns['mul'](angle.cos(),ns['direct']),ns['mul'](angle.sin(),ns['normal']));u=ns['mul'](-1,ns['perp'](direction));slope=ns['dot'](u,ns['D']);assert slope<0
 selected=[x for x in doc['dihedral_witnesses'] if x['edge']=='AH-cut'];assert len(selected)==2 and all(arb(x['violation'])>0 for x in selected)
 doc.update(pivot_exclusion='all p >= 83331/100000',cut_offset_derivative=str(slope),input_sha256=sha256(cfgpath.read_bytes()).hexdigest(),replay_wrapper_sha256=sha256(Path(__file__).read_bytes()).hexdigest());a.output.write_text(json.dumps(doc,indent=2)+'\n');raw.unlink();print('pivot >=',cfg['pivot'],'excluded; slope',slope)
if __name__=='__main__':main()
