#!/usr/bin/env python3
"""Certify EH root existence by four ordered branches of a two-line ear."""
import argparse
from fractions import Fraction
from itertools import product
import json
from pathlib import Path

import numpy as np

import audit_optimized_eh_chord as g
from certify_eh_joint_switch import add_joint_query, certify
from certify_eh_joint_parameters import gain_endpoint, specialize
from certify_eh_polyline_obstruction import clip_geometry
from certify_eh_relocation_root_condition import phi_objective
from certify_eh_sprague_relocation import arc_data, voronoi_signs


def run(parts, right, max_depth, receipt=None):
    sigma = Fraction('1.5497')
    constructor, model_hash = specialize(sigma)
    normals, offsets, gain, _, _ = clip_geometry(parts)
    outer_depth = parts[0]
    if receipt:
        assert receipt['parameters'] == parts and Fraction(receipt['right_degrees']) == right
        assert receipt['specialized_model_sha256'] == model_hash

    branch_certs = {}
    bounds = []
    for pair in product(range(2), repeat=2):
        program = constructor(outer_depth)
        coordinates = np.eye(program.dimension+1, dtype=object)
        for i, edge in enumerate(pair):
            candidate = coordinates[2*i:2*i+2]
            program.linear(-np.asarray(normals[edge], dtype=object)@candidate, -offsets[edge])
        key = ''.join(map(str, pair))
        objective = phi_objective(program, right, 1)
        if receipt:
            witness = receipt['ordered_branches'][key]
            for i, direction in witness['adaptive_cuts']:
                program.cut(i, direction, record=True)
            weights = witness['dual']
        else:
            weights = program.optimize(objective)
        estimate = program.validate(objective, weights)
        assert estimate < g.arb('1.5'), (key, str(estimate-g.arb('1.5')))
        branch_certs[key] = {'F_right_upper': str(estimate-g.arb('1.5')), 'dual': weights,
                             'adaptive_cuts': program.adaptive_cuts}
        bounds.append(estimate)

    # This valid inequality follows for the same body in any of the four branches.
    # All three coefficients of the support sum are positive, so it also holds
    # for the arbitrary support witnesses used in the enclosing-ear model.
    enclosing_program = constructor(outer_depth)
    add_joint_query(enclosing_program)
    enclosing_program.linear(phi_objective(enclosing_program, right, 1), max(bounds))
    saved = receipt['enclosing_pose_certificate'] if receipt else None
    if saved:
        for i, direction in saved['adaptive_cuts']:
            enclosing_program.cut(i, direction, record=True)
    placement = certify(enclosing_program, right, max_depth, saved, eh_relocation=True,
                        joint_eh=True, revisit_ah=True)

    _, _, c3, b3, _, _, _, _, j = g.setup()
    v, w = g.minus(j, c3), g.minus(j, b3)
    assert g.cross(v, w) > 0 and g.dot(v, w) > 0
    area = [str(gain_endpoint(b)) for b in (gain.lower(), gain.upper())]
    composition_summary = {**placement['summary'], 'sigma_degrees': str(sigma),
                           'delta_is_enclosing_placement_depth': True,
                           'cut_gain_area_interval_exact': area}
    return {'parameters': parts, 'right_degrees': str(right), 'sigma_degrees': str(sigma),
            'specialized_model_sha256': model_hash,
            'cut_gain_area_interval_exact': area, 'summary': composition_summary,
            'sprague_endpoint_signs': voronoi_signs(arc_data()),
            'ordered_branches': branch_certs, 'enclosing_pose_certificate': placement}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('first_depth')
    parser.add_argument('second_depth')
    parser.add_argument('hinge_fraction')
    parser.add_argument('output', type=Path)
    parser.add_argument('--right', default='36.01')
    parser.add_argument('--max-depth', type=int, default=10)
    parser.add_argument('--replay', type=Path)
    args = parser.parse_args()
    frozen = json.loads(args.replay.read_bytes()) if args.replay else None
    data = run([args.first_depth, args.second_depth, args.hinge_fraction], Fraction(args.right), args.max_depth, frozen)
    args.output.write_text(json.dumps(data, separators=(',', ':'))+'\n')
