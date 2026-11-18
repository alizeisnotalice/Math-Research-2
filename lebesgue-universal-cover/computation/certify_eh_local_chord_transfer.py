"""Combine the R176 three-facet ear with an AH certificate using its true EH cap rectangle."""
import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path

import certify_eh_third_facet as base

g=base.g


def verify(archive,ah_file,eh_file):
    ah=json.loads(ah_file.read_bytes());eh=json.loads(eh_file.read_bytes())
    r126=archive/'eh-two-line-cover-certificate.json';old=json.loads(r126.read_bytes())
    checksum=sha256(r126.read_bytes()).hexdigest()
    assert eh['R126_input_sha256']==ah['source_hashes'][r126.name]==checksum
    assert eh['depth']==base.DEPTH and F(eh['hinge_fraction_exact'])==base.PIVOT
    _,model_hash=base.eh.specialize(F('1.5497')); assert eh['model_sha256']==model_hash
    assert eh['weights_sha256']==sha256((archive/'eh-third-facet-weights.json').read_bytes()).hexdigest()
    assert len(eh['ordered_F_upper'])==9 and all(g.arb(x)<base.arbq(base.F_LIMIT) for x in eh['ordered_F_upper'].values())
    assert eh['EH_placement']['angle_edge_cells']==54 and eh['EH_placement']['validated_duals']==63
    _,_,points,u,gain=base.geometry(old);v,_,_,t=points
    assert v[0]<t[0] and t[1]<v[1]
    n_u,n_v=base.unit_outer(u,t),base.unit_outer(v,t)
    # A point violating the deeper chord lies inside the closed two-arc cap
    # of the true pre-EH host. Its coordinate box has T as a corner, so it
    # suffices to pay the norm change on this box, including both arcs.
    change=g.norm(g.minus(n_u,n_v));diameter=g.norm(g.minus(v,t));cost=change*diameter
    placement_rows = [*ah['bins'], *ah.get('root_conditioned_classes', [])]
    margin=min(g.arb(row['placement_margin_lower']).lower() for row in placement_rows);assert margin>cost
    corners=[(a,b) for a in (v[0],t[0]) for b in (v[1],t[1])]
    clear=base.prev.ah_clearance(ah,corners)
    gl,gu=map(F,eh['EH_gain_over_R126_exact_interval']);assert (gl,gu)==tuple(map(base.eh.gain_endpoint,(gain.lower(),gain.upper())))
    al,au=map(F,ah['area_exact_interval']);bounds=al-gu,au-gl
    return {'depth':base.DEPTH,'pivot':str(base.PIVOT),'AH_input_sha256':sha256(ah_file.read_bytes()).hexdigest(),
            'EH_input_sha256':sha256(eh_file.read_bytes()).hexdigest(),'old_AH_placement_minimum_lower':str(margin),
            'rectangle_diagonal_upper':str(diameter.upper()),'unit_normal_change_upper':str(change.upper()),
            'local_chord_transfer_cost_upper':str(cost.upper()),'AH_after_EH_transfer_margin_lower':str((margin-cost).lower()),
            'AH_rectangle_clearance_lower':clear,'EH_gain_over_R126_exact_interval':[str(gl),str(gu)],
            'area_exact_interval':list(map(str,bounds)),'area_decimal_interval':base.prev.format_area(*bounds)}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('archive',type=Path);p.add_argument('AH',type=Path);p.add_argument('EH',type=Path);p.add_argument('output',type=Path)
    a=p.parse_args();a.output.write_text(json.dumps(verify(a.archive,a.AH,a.EH),indent=2)+'\n')
