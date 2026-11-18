"""Replay the eleven-facet EH certificate with direct AH root tests of its deepest chord."""
import argparse
from fractions import Fraction as F
from hashlib import sha256
import inspect,json,sys
from pathlib import Path

from flint import arb,ctx


def run(archive,ah_path,parameters,weights,output):
    from certify_eh_upper_chain import run as original
    from certify_ah_three_line_cut import format_area
    from certify_ah_direct_eh_row import add_row

    expected=archive/'eh-upper-chain-parameters.json'
    assert parameters.read_bytes()==expected.read_bytes()
    ah=json.loads(ah_path.read_bytes())
    ctx.prec=128
    observed=add_row({'rows':[]},archive)
    assert observed==ah['direct_eh_chord']
    assert ah['complete'] and ah['family']['direct_eh_chord'] is True
    assert observed['number_of_EH_facets']==11

    src=inspect.getsource(original)
    first=src.index('    m=placement_lower(ah)')
    last=src.index('    clear=b.prev.ah_clearance',first)
    src=src[:first]+'''    slack=placement_lower(ah)
    direct=[]
    assert ah['direct_eh_chord']['parameters_sha256']==sha256(param_file.read_bytes()).hexdigest()
    assert ah['direct_eh_chord']['R126_sha256']==sha256(oldpath.read_bytes()).hexdigest()
    assert len(geo['active'])==ah['direct_eh_chord']['number_of_EH_facets']==11
    assert slack>0
'''+src[last:]
    env=dict(original.__globals__)
    exec(compile(src,'eleven_EH_with_direct_AH_root_row','exec'),env)
    result=env['run'](archive,ah_path,parameters,json.loads(weights.read_bytes()),False,None)
    result['weights_sha256']=sha256(weights.read_bytes()).hexdigest()
    result['direct_EH_row_sha256']=sha256(src.encode()).hexdigest()
    result['direct_EH_chord_checked_in_all_AH_placements']=True

    eh_low,eh_high=map(F,result['EH_gain_over_R126_exact_interval'])
    sl,su=map(F,ah['strong_member_area_exact'])
    wl,wu=map(F,ah['weak_member_area_exact'])
    assert [sl-eh_high,su-eh_low]==list(map(F,result['area_exact_interval']))
    assert sl<=su<wl<=wu
    uniform=min(arb(ah['uniform_cut_line_margin_lower']).lower(),
                arb(result['AH_after_transfer_lower']).lower())
    assert uniform>0
    result['all_members_area_exact_interval']=list(map(str,(sl-eh_high,wu-eh_low)))
    result['all_members_area_decimal_interval']=format_area(sl-eh_high,wu-eh_low)
    result['weak_member_area_decimal_interval']=format_area(wl-eh_high,wu-eh_low)
    result['uniform_direct_AH_and_EH_slack_lower']=str(uniform.lower())
    output.write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('archive','AH','parameters','weights','out'):
        p.add_argument(k,type=Path)
    a=p.parse_args()
    archive=a.archive.resolve()
    sys.path.insert(0,str(archive))
    run(archive,a.AH,a.parameters,a.weights,a.out)
