"""Replay the first EH placement and certify transfer to the real 1.54 degree AH chart."""
from argparse import ArgumentParser
from fractions import Fraction as F
from hashlib import sha256
import json
import os
from pathlib import Path
import subprocess
import sys

from flint import arb,ctx

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
from certify_three_sigma_area_budget import fixed_ah,pack,subtract,setup,ah_clearance

SIGMA='1.54'
PIVOT=F('.429')
TURN=F('.0062')
D_CHART=F('.000555')
DEPTH='0.000550'
RIGHT='36.01'


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def load(path):
    return json.loads(path.read_bytes())


def exact(q):
    q=F(q)
    return arb(q.numerator)/q.denominator


def replay_base(comp,work):
    env={**os.environ,'LEBESGUE_SIGMA_DEGREES':'1.54,1.54'}
    sources=('certify_ah_shortest.py','certify_ah_contact_order.py','audit_eh_geometry_final.py',
             'certify_splice_connectors.py','certify_cover_relocations.py')
    fingerprints={}
    for source in sources:
        out=subprocess.check_output([sys.executable,'-B',str(comp/source)],env=env)
        data=json.loads(out)
        assert data['sigma_degrees_interval']==['77/50','77/50']
        assert data.get('passed',data.get('exact_maps_passed')) is True
        fingerprints[source]=sha256(out).hexdigest()
    raw=subprocess.check_output([sys.executable,'-B',str(comp/'certify_gibbs_area_refined.py'),SIGMA,'8','32'])
    area_path=work/'old-area.json'
    area_path.write_bytes(raw)
    area=load(area_path)
    assert F(area['sigma_degrees'])==F(SIGMA) and area['passed']
    interval=tuple(map(F,area['area_interval']))
    return area_path,interval,fingerprints


def run(a):
    comp=a.comp.resolve()
    a.work.mkdir(parents=True,exist_ok=True)
    sys.path.insert(0,str(comp))
    import audit_optimized_eh_chord as g
    from certify_sigma_ah_full_chart import geometry
    from certify_sigma_reachable_roots import reach
    from certify_eh_joint_parameters import specialize

    basefile,base_area,base_geometry=replay_base(comp,a.work)
    data=load(a.eh)
    before=data['summary']
    assert F(before['sigma_degrees'])==F(SIGMA) and F(before['delta'])==F(DEPTH)
    reply=a.work/'eh-replayed.json'
    command=[sys.executable,'-B',str(comp/'certify_eh_joint_switch.py'),
             SIGMA,DEPTH,RIGHT,str(reply),'--max-depth','9',
             '--eh-relocation','--joint-eh','--revisit-ah',
             '--old-area',str(basefile),'--replay',str(a.eh)]
    subprocess.run(command,check=True,stdout=sys.stderr)
    fresh=load(reply)
    after=fresh['summary']
    assert after==before
    assert arb(after['F_right_upper'])<0 and arb(after['minimum_slack_lower_bound'])>0
    assert after['test_edge_count']==4 and F(after['root_angle_degrees'][1])==F(RIGHT)
    assert after['AH_relocation_after_root']

    ahr=load(a.ah)
    cfg=load(a.config)
    assert ahr['input_sha256']==digest(a.config) and ahr['parameters']==cfg
    assert ahr['source_sha256']==digest(comp/'certify_sigma_reachable_roots.py')
    assert F(cfg['sigma_degrees'])==F(SIGMA)
    assert list(map(F,cfg['line']))==[PIVOT,TURN]
    far,negative,eps=map(F,cfg['physical_parameters'])
    ns,src,_=geometry(comp,[(PIVOT,TURN)],sigma=SIGMA,far=far,negative=negative,eps=eps)
    assert sha256(src.encode()).hexdigest()==ahr['physical_geometry_sha256']
    assert F(ns['eh_depth'])==D_CHART
    active=[j for j,E in enumerate(ns['deficit_by_bin']) if not E.is_zero()]
    assert active==[r[0] for r in cfg['bins']]==[r['reference_bin'] for r in ahr['rows']]==list(range(70))
    assert ahr['all_displayed_bins_passed']
    root_signs=[];margin=[]
    for row in ahr['rows']:
        cand=[list(z) for z in reach(ns,row['reference_bin'])]
        assert row['passed'] and row['reachable_first_hit_subcell_ranges']==cand
        for trial in row['tested']:
            root_signs.append(arb(trial['root_margin_lower']))
            margin.append(arb(trial['placement_margin_lower']))
    assert min(root_signs)>0 and min(margin)>0
    assert F(DEPTH)<D_CHART

    # Check the two EH chords and the actual retained side at the chart depth.
    ctx.prec=160
    _,_=specialize(F(SIGMA))
    V=setup(SIGMA,g)
    U,C3,T=[V[k] for k in ('U','C3','T')]
    radial=g.minus(U,C3)
    def rotated(depth):
        return g.plus(C3,g.scale(1,g.plus(g.scale(exact(depth).cos(),radial),
                                          g.scale(-exact(depth).sin(),g.rot(radial)))))
    chart_u=rotated(D_CHART)
    shallower=rotated(F(DEPTH))
    u,v=g.minus(chart_u,C3),g.minus(shallower,C3)
    assert g.cross(u,v)>0 and g.dot(u,v)>0 and chart_u[1]>shallower[1]>T[1]
    cut=g.minus(T,chart_u)
    chart_out=g.scale(-1/g.norm(cut),g.rot(cut))
    chart_offset=g.dot(chart_out,T)
    (label,n,b),=[z for z in ns['rows'] if z[0].startswith('U1-Te')]
    assert all(x.overlaps(y) for x,y in zip(n,chart_out,strict=True)) and b.overlaps(chart_offset)
    cut_sign=-g.cross(cut,g.minus(shallower,chart_u))
    assert cut_sign>0
    clearance=ah_clearance(V,shallower,T,g,PIVOT,TURN)
    ah_gain,w,slope=fixed_ah(comp,SIGMA,PIVOT,TURN,far)
    assert 0<w[0]<w[1]<far
    # The paired AH edge is transverse at the far endpoint; old three
    # shallower AH lines are pointwise weaker in the only active w range.
    for p,d in ns['cuts'][:-1]:
        assert p<PIVOT and d<TURN
    assert ns['cuts'][-1]==(PIVOT,TURN)

    eh_gain=tuple(map(F,after['cut_gain_area_interval_exact']))
    total=subtract(subtract(base_area,ah_gain),eh_gain)
    r201=comp/'fan-0062-eleven-eh.json'
    prev=tuple(map(F,load(r201)['area_exact_interval']))
    saving=subtract(prev,total)
    assert saving[0]>0
    warning=after['active_edge']
    return dict(sigma=SIGMA,pivot=str(PIVOT),AH_turn=str(TURN),EH_depth=DEPTH,
                EH_chart_depth=str(D_CHART),EH_scope=after['root_angle_degrees'],
                EH_root_upper=after['F_right_upper'],
                EH_placement_margin=after['minimum_slack_lower_bound'],
                EH_weakest_dir=warning,EH_cells=after['angle_edge_cells'],
                EH_duals=after['validated_duals'],charts=len(active),
                old_geometry_replays=base_geometry,
                AH_root_margin=pack(min(root_signs)),AH_placement_margin=pack(min(margin)),
                chart_to_single_EH_nested_determinant=pack(cut_sign),
                AH_EH_clearance=pack(clearance),AH_curve_crossing=pack(w),
                AH_curve_crossing_deriv_lower=str(slope),
                base_area=pack(base_area),AH_cut=pack(ah_gain),EH_cut=pack(eh_gain),
                combined_area=pack(total),saving_from_R201=pack(saving),
                eh_receipt_sha256=digest(a.eh),ah_receipt_sha256=digest(a.ah),
                ah_input_sha256=digest(a.config),source_sha256=digest(Path(__file__)),
                later_reference_sha256=digest(r201),
                replayed_eh_sha256=digest(reply),
                precision_bits=ctx.prec)


if __name__=='__main__':
    assert __debug__
    parser=ArgumentParser(description=__doc__)
    for key in ('comp','eh','ah','config','work'):
        parser.add_argument(key,type=Path)
    print(json.dumps(run(parser.parse_args()),separators=(',',':')))
