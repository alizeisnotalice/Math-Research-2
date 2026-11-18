"""Verify a common high root with early reflection for the entire bin 8 of a two-chord cap."""
import argparse
from fractions import Fraction as F
from hashlib import sha256
import inspect,json,sys
from pathlib import Path
from flint import arb
from certify_ah_coupled_brackets import scan
from certify_ah_margin_envelope import replace_line
from certify_ah_wide_contact_chart import geometry


def run(archive,inp,out):
    cfg=json.loads(inp.read_bytes());lines=[tuple(map(F,r)) for r in cfg['lines']];ns,source=geometry(archive,lines)
    code=inspect.getsource(scan);code=replace_line(code,'(300,1,1)','(120,-1,1)');sc=dict(scan.__globals__);exec(compile(code,'early_reflected_high_root','exec'),sc);result=sc['scan'](ns,8,0,75,*map(F,cfg['angles_degrees']),32)
    assert arb(result['root_slack_lower'])>0 and arb(result['placement_slack_lower'])>0
    out.write_text(json.dumps({'input':cfg,'pose':{'turn_degrees':120,'orientation_sign':-1,'root_sign_at_left':1},'parameter_cells':75,'angle_cells':2400,'geometry_sha256':sha256(source.encode()).hexdigest(),'derived_scan_sha256':sha256(code.encode()).hexdigest(),'input_sha256':sha256(inp.read_bytes()).hexdigest(),'script_sha256':sha256(Path(__file__).read_bytes()).hexdigest(),'common_contact_kernel_sha256':sha256(Path(inspect.getfile(scan)).read_bytes()).hexdigest(),'python':sys.version,'flint':ns['flint_version'],'result':result},indent=2)+'\n')

if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('archive',type=Path);ap.add_argument('input',type=Path);ap.add_argument('output',type=Path);opt=ap.parse_args();sys.path.insert(0,str(opt.archive.resolve()));run(opt.archive.resolve(),opt.input,opt.output)
