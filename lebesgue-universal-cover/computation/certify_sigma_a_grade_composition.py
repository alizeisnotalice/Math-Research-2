"""Bind real double-hit placement, the two initial EH branches, and signed true Green areas at 1.54 degrees."""
import argparse,json,sys
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
from flint import arb,ctx
from certify_ah_three_line_cut import format_area


def run(ah,eh,records,area_path,out):
    ctx.prec=160;hashes={}
    def read(p):
        raw=p.read_bytes();hashes[p.name]=sha256(raw).hexdigest();return json.loads(raw)
    joined=read(ah/'sigma-154-a-grade-complete-ah.json');single=read(eh/'sigma154-ah429-eh550-composition.json');eleven=read(eh/'sigma154-ah429-eh-chain997-composition.json');area=read(area_path);green=read(ah/'sigma-154-full-budget-certificate.json');hist=read(records/'flexible-006-eleven-eh.json');latest=read(records/'fan-0070-final-eleven-eh.json')
    assert joined['number_reference_bins']==70 and joined['reachable_reference_subcells']==5181 and arb(joined['root_margin_lower'])>0
    assert [F(x) for x in (area['sigma'],area['pivot'],area['turn'])]==list(map(F,['1.54','.429','.009']))
    assert area['single_sha256']==hashes['sigma154-ah429-eh550-composition.json'] and area['eleven_sha256']==hashes['sigma154-ah429-eh-chain997-composition.json']
    assert F(area['AH_curve_crossing']['exact'][1])<F(area['AH_far'])==F('.273') and F(area['AH_crossing_derivative_lower'])>0
    assert arb(single['EH_root_upper'])<0 and arb(single['EH_placement_margin'])>0 and F(single['chart_to_single_EH_nested_determinant']['exact'][0])>0
    assert F(eleven['EH_right_critical_upper']['exact'][1])<0 and arb(eleven['EH_pose_margin'])>0 and sorted(eleven['EH_11_active_facets'])==list(range(11))
    transfer=F(eleven['transfer_cost']['exact'][1]);cost=arb(transfer.numerator)/transfer.denominator;room=arb(joined['placement_margin_lower'])-cost;assert room>0
    b0,b1=map(F,green['gibbs_base']['area_interval']);arow=next(r for r in green['rows'] if list(map(F,r['AH_chord']))==list(map(F,['.429','.009'])));h0,h1=map(F,arow['AH_gain_exact_interval']);ha,hb=map(F,area['AH_cut']['exact']);assert max(h0,ha)<min(h1,hb);h0=max(ha,h0);h1=min(hb,h1)
    rows=[];old0,old1=map(F,hist['area_exact_interval']);prev0,prev1=map(F,latest['area_exact_interval'])
    for entry in area['cases']:
        assert F(entry['AH_EH_clearance']['exact'][0])>0
        e0,e1=map(F,entry['EH_cut']['exact']);L,U=b0-h1-e1,b1-h0-e0;coarse0,coarse1=map(F,entry['resulting_area']['exact']);assert coarse0<L<U<coarse1 and old0-U>F(1,10**6) and prev0>U
        rows.append({'EH':entry['label'],'area_exact_interval':list(map(str,(L,U))),'area_decimal_interval':format_area(L,U),'saving_from_round_199':format_area(old0-U,old1-L),'saving_from_round_203':format_area(prev0-U,prev1-L),'beyond_A_margin_lower':str(old0-U-F(1,10**6))})
    obj={'sigma_degrees':'1.54','AH_pivot':'.429','AH_angle_radians':'.009','initial_hexagon_sign_and_placement':{'single_root':single['EH_root_upper'],'single_pose':single['EH_placement_margin'],'eleven_root':eleven['EH_right_critical_upper'],'eleven_pose':eleven['EH_pose_margin']},'actual_AH_root_lower':joined['root_margin_lower'],'actual_AH_all_edges_lower':joined['placement_margin_lower'],'deducted_transfer_upper':str(transfer),'retained_all_edge_margin':str(room.lower()),'covers':rows,'input_sha256':hashes,'source_sha256':sha256(Path(__file__).read_bytes()).hexdigest(),'python':sys.version};out.write_text(json.dumps(obj,indent=2)+'\n');print(rows)


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('ah',type=Path);ap.add_argument('eh',type=Path);ap.add_argument('records',type=Path);ap.add_argument('area',type=Path);ap.add_argument('out',type=Path);o=ap.parse_args();run(o.ah,o.eh,o.records,o.area,o.out)
