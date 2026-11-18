"""Replay both AH root branches and the exact area of a two-line cap."""
import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
import sys

from flint import arb, ctx


def run(archive, inp, out):
    sys.path.insert(0,str(archive))
    from certify_ah_wider_dual_chord import prepare, area_and_witnesses
    from wider_chart import bracket
    from certify_ah_early_reflection import run as replay_early
    from certify_ah_three_line_cut import format_area

    config=json.loads(inp.read_bytes())
    lines=config['lines']
    ar,source,exact_lines=prepare(archive,lines)
    lo,hi,shape,witnesses=area_and_witnesses(archive,ar,exact_lines)
    ctx.prec=128
    required={j for j,E in enumerate(ar['deficit_by_bin']) if not E.is_zero()}
    assert 8 in required
    general,program=bracket(ar,required-{8})
    assert {r['bin'] for r in general}==required-{8}
    early_path=out.with_name('early-reflection.json')
    replay_early(archive,inp,early_path)
    early=json.loads(early_path.read_bytes())
    block=early['result']
    geom_sha=sha256(source.encode()).hexdigest()
    assert early['geometry_sha256']==geom_sha
    assert block['reference_bin']==8
    assert list(map(F,block['closed_w_endpoints']))==[F('0.08'),F('0.09')]
    assert block['subcell_halfopen_indices']==[0,early['parameter_cells']]
    assert arb(block['root_slack_lower'])>0 and arb(block['placement_slack_lower'])>0
    entry={'bin':8,'parameter_cells':early['parameter_cells'],'angle_cells':early['angle_cells'],
           'placement_margin_lower':block['placement_slack_lower'],'root_margin_lower':block['root_slack_lower']}
    all_bins=sorted([*general,entry],key=lambda r:r['bin'])
    value={'complete':True,'scope':'all constant-width-one inputs by the two AH root branches',
           'lines':lines,'bins':all_bins,'shape':shape,'R150_R156_witnesses_at_identity':witnesses,
           'exceptional_bin_proof':early,'geometry_sha256':geom_sha,
           'ordinary_bin_verifier_sha256':sha256(program.encode()).hexdigest(),
           'source_hashes':{name:sha256((archive/name).read_bytes()).hexdigest() for name in ('eh-two-line-cover-certificate.json','gibbs-1p5497-refined-18-128.json','certify_ah_combined_area.py')},
           'R167_area_sha256':sha256((archive/'ah-three-line-cut-certificate.json').read_bytes()).hexdigest(),
           'area_exact_interval':[str(lo),str(hi)],'area_decimal_interval':format_area(lo,hi),
           'total_parameter_cells':sum(r['parameter_cells'] for r in all_bins),
           'total_angle_cells':sum(r['angle_cells'] for r in all_bins),
           'python':sys.version,'flint':ar['flint_version']}
    out.write_text(json.dumps(value,indent=2)+'\n')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ('archive','input','output'):
        parser.add_argument(name,type=Path)
    args=parser.parse_args()
    run(args.archive,args.input,args.output)
