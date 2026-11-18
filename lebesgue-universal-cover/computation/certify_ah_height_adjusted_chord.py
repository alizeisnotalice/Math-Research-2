"""Physical paired-height debit, full placements and Green area for a single deep AH chord."""
import argparse,inspect,json,sys
from fractions import Fraction as F
from pathlib import Path
from hashlib import sha256
from flint import arb
from certify_ah_margin_envelope import replace_line
from certify_ah_paired_height import geometry
from certify_ah_projected_jets import bracket
from certify_ah_root_conditioned_support import root_conditioned
import certify_ah_deep_single_line as source


def run(archive,inp,out):
    cfg=json.loads(inp.read_bytes()); (p,d),=[tuple(map(F,r)) for r in cfg['lines']];far,negative,eps=map(F,cfg['physical_parameters']);ns,full,_=geometry(archive,[(p,d)],far=far,negative=negative,eps=eps)
    explicit={j for low,high,_,_ in cfg['common_brackets'] for j in range(low,high)}
    remaining=[j for j,e in enumerate(ns['deficit_by_bin']) if not e.is_zero() and j not in explicit]
    regular,verify_source=bracket(ns,remaining); common=[]
    for low,high,l,r in cfg['common_brackets']:
        for j in range(low,high):
            record=root_conditioned(ns,j,F(l),F(r),16);assert arb(record['root_margin_lower'])>0 and arb(record['placement_margin_lower'])>0,record
            common.append(record)
    assert len(common)==len(explicit) and {r['reference_bin'] for r in common}.isdisjoint(row['bin'] for row in regular)
    # Reuse the single-chord area and guard checks; the first phase just completed for this exact chart.
    def reconstruct(path,extra):
        assert path==archive and extra==[(p,d)];return ns,full
    def checked(_ns,_bins):
        assert _ns is ns and _bins is None;return regular,verify_source
    green_src=inspect.getsource(source.green);green_src=replace_line(green_src,'F("0.07")',f'F("{far}")');gn=dict(source.green.__globals__);exec(compile(green_src,'height_adjusted_green','exec'),gn)
    driver=inspect.getsource(source.run);driver=replace_line(driver,'json.loads(inp.read_bytes())',"json.loads(inp.read_bytes())['lines']") if driver.count('json.loads(inp.read_bytes())')==1 else driver.replace('json.loads(inp.read_bytes())',"json.loads(inp.read_bytes())['lines']")
    driver=replace_line(driver,"0<d<=F('.0006')","0<d<F('.004')");driver=replace_line(driver,"hi<F('.844093390960250720226452027837524452015')","hi<F('.844093233635451252504326017398854271746')")
    stages={**source.run.__globals__,'make_geometry':reconstruct,'fine_verify':checked,'green':gn['green'],'__file__':__file__};exec(compile(driver,'adjusted_chord_driver','exec'),stages);stages['run'](archive,inp,out)
    data=json.loads(out.read_bytes());data['physical_parameters']=cfg['physical_parameters'];data['root_conditioned_classes']=common;data['debit_kernel_sha256']=sha256(Path(__file__).with_name('certify_ah_paired_height.py').read_bytes()).hexdigest();out.write_text(json.dumps(data,indent=2)+'\n')

if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('archive',type=Path);ap.add_argument('input',type=Path);ap.add_argument('output',type=Path);opt=ap.parse_args();sys.path.insert(0,str(opt.archive.resolve()));run(opt.archive.resolve(),opt.input,opt.output)
